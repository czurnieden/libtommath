#include "tommath_private.h"
#ifdef BN_MP_GET_STR_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */




/*
   Old prototype of a fast number-conversion for LibTomMath
   Needs fast division to get all of its potential, of course,
   but works without, too.
*/

#include <string.h>
#include <stdint.h>

static int s_ilog2(int value)
{
   int r = 0;
   while ((value >>= 1) != 0) {
      r++;
   }
   return r;
}

static void *s_mp_memset(void *s, int c, size_t n)
{
   char *_s = (char *) s;
   while (n-- != 0u) {
      *_s = (char) c;
      _s++;
   }
   return s;
}

static int s_mp_strlen(const char *s)
{
   const char *p;
   p = s;
   while (*p != '\0') {
      p++;
   }
   return (int)(p - s);
}

char *s_mp_strncat(char *s1, const char *s2, size_t n)
{
   char *s = s1;

   while (*s1 != '\0') {
      s1++;
   }
   while (n-- != 0) {
      *s1 = *s2;
      if (*s2 == '\0') {
         break;
      }
      s1++;
      s2++;
      if (n == 0) {
         *s1 = '\0';
      }
   }
   return s;
}


/* Size is in bits (75 bytes for base 2 is comfortably small enough)*/
#define SCHOENHAGE_CONVERSION_CUT_BASE  600
/* There is a possibility of a bit of tune-ability here, I may note. */
#define SCHOENHAGE_CONVERSION_CUT  600

/* based on SCHOENHAGE_CONVERSION_CUT=600 */
/*
   Computed with (Pari/GP)
      len(cut, radix) = ceil(log(2^(cut+1))/log(radix)) + 1
      for(i = 2, 64, k = len(600,i);printf(k ", "))

    Although the table is hardcoded for 600 bit a change to SCHOENHAGE_CONVERSION_CUT
    is directly proportional to the sizes in the table, e.g.: a raise of
    the SCHOENHAGE_CONVERSION_CUT value by 10% gives rise to the values in the table by
    10%, too. To add some angst-allowance for the unavoidable rounding errors is
    recommended nevertheless.
*/

static const size_t buf_size[] = {
   0, 0, 602, 381, 302, 260, 234, 216, 202, 191, 182, 175, 169, 164, 159, 155, 152,
   149, 146, 143, 141, 138, 136, 134, 133, 131, 129, 128, 127, 125, 124, 123,
   122, 121, 120, 119, 118, 117, 116, 115, 114, 114, 113, 112, 112, 111, 110,
   110, 109, 109, 108, 107, 107, 106, 106, 105, 105, 105, 104, 104, 103, 103,
   102, 102, 102
};


/* floor(log_2(x)) */
static const int log_table[65] = {
   0, 0, 1, 1, 2, 2, 2, 2,
   3, 3, 3, 3, 3, 3, 3, 3,
   4, 4, 4, 4, 4, 4, 4, 4,
   4, 4, 4, 4, 4, 4, 4, 4,
   4, 5, 5, 5, 5, 5, 5, 5,
   5, 5, 5, 5, 5, 5, 5, 5,
   5, 5, 5, 5, 5, 5, 5, 5,
   5, 5, 5, 5, 5, 5, 5, 5,
   6
};

static mp_err s_mp_get_str_intern(mp_int *a, char *string, int digits, int base,
                                  size_t size, mp_int *s_schoenhagecache)
{
   int b, n, i;
   int ed;
   mp_err err;
   mp_int q, r;
   char *_s, buf[size + 1];

   /* Use naive method for the leaves */
   if ((a->used * MP_DIGIT_BIT) < SCHOENHAGE_CONVERSION_CUT) {
      if ((err = mp_to_radix(a, buf, SIZE_MAX, base)) != MP_OKAY) {
         return err;
      }
      /*
          I think that there is a slight chance to put those leading zeros
          in their places in a more elegant way ;-)
       */
      /* Can be replaced with the value form the new mp_to_radix */
      i = s_mp_strlen(buf);
      if ((i < digits) && (*string != '\0')) {
         _s = string;
         while (*string != '\0') {
            string++;
         }
         for (; i < digits; i++) {
            *string = (char) '0';
            string++;
         }
         string = _s;
      }
      string = s_mp_strncat(string, buf, SIZE_MAX);
      return MP_OKAY;
   }

   b = mp_count_bits(a);
   /* Compute the exponent for the base necessary */
   n = s_ilog2(b / (2 * log_table[base])) - 1;
   if (n >= (int)(sizeof(int) * CHAR_BIT) - 1) {
      return MP_VAL;
   }

   if ((err = mp_init_multi(&q, &r, NULL)) != MP_OKAY) {
      return err;
   }
   /*
       Divide a chunk of the input by the correct basepower from the cache.
       This is the point where a fast division algorithm would come handy.
    */
   if ((err = mp_div(a, &(s_schoenhagecache[n]), &q, &r)) != MP_OKAY)                goto LBL_ERR;
   ed = 1 << n;
   /* Rinse and repeat with quotient and remainder */
   if ((err = s_mp_get_str_intern(&q, string, digits - ed,
                                  base, size, s_schoenhagecache)) != MP_OKAY)    goto LBL_ERR;
   if ((err = s_mp_get_str_intern(&r, string, ed,
                                  base, size, s_schoenhagecache)) != MP_OKAY)             goto LBL_ERR;

LBL_ERR:
   mp_clear_multi(&q, &r, NULL);
   return MP_OKAY;
}

/* The buffer "string" must be allocated in advance */
/* TODO: do it internally when the new mp_to_radix comes out */
mp_err mp_get_str(const mp_int *a, char *string, size_t maxlen, size_t *written, int base)
{
   mp_sign sign;
   mp_err err = MP_OKAY;
   mp_int a_bis;

   int s_schoenhagecache_len = 0, i, n, b;
   size_t buffer_size, wrote = 0u;

   mp_int s_schoenhagecache[(sizeof(int) * CHAR_BIT) * sizeof(mp_int)];

   /* check range of the maxlen, radix */
   if ((maxlen < 2u) || (base < 2) || (base > 64)) {
      return MP_VAL;
   }

   if (string == NULL) {
      return MP_VAL;
   }

   /* quick out if its zero */
   if (MP_IS_ZERO(a)) {
      *string++ = '0';
      *string = '\0';
      if (written != NULL) {
         *written = 2;
      }
      return MP_OKAY;
   }

   b = mp_count_bits(a);
   if (b <= (SCHOENHAGE_CONVERSION_CUT)) {
      /* TODO: adapt to new mp_to_radix */
      if ((err = mp_to_radix(a, string, SIZE_MAX, base)) != MP_OKAY) {
         return err;
      }
      if (written != NULL) {
         /* TODO: remove after new mp_to_radix is included */
         wrote = (size_t)s_mp_strlen(string) + 1;
         *written = wrote;
      }
   }
   /*
      TODO: To avoid unnecessary complications we check the output-buffer size here.
            For now!
            Caveat: as mp_radix_size_overestimate does indeed overestimate

               if (maxlen < size) {
                  return MP_VAL;
               }

            will not work in that case!

   if ((err = mp_radix_size(a, base, &size) ) != MP_OKAY){
      return err;
   }
   if (maxlen < (size_t)size) {
      return MP_VAL;
   }
   */
   /* Compute size of buffer */
   buffer_size = (buf_size[base] * SCHOENHAGE_CONVERSION_CUT) / SCHOENHAGE_CONVERSION_CUT_BASE;

   /* Compute the max exponent for the base necessary */
   n = s_ilog2(b / (2 * log_table[base])) - 1;
   if (n >= (int)(sizeof(int) * CHAR_BIT) - 1) {
      return MP_VAL;
   }
   /*
       Fill the cache with squares of squares of "base". Example with base 10 (ten)

       s_schoenhagecache[0] = 10
       s_schoenhagecache[1] = s_schoenhagecache[0]^2 = 100
       s_schoenhagecache[2] = s_schoenhagecache[1]^2 = 10000
       ...

       TODO: we don't need the entries of size < (SCHOENHAGE_CONVERSION_CUT)
             (in their respective bases, of course)

       TODO: The entries get used multiple times, so it is a good idea to replace the entries
             with their approximate reciprocals and use Barrett division.

    */
   if ((err = mp_init_set(&(s_schoenhagecache[0]), (mp_digit)(base))) != MP_OKAY) goto LBL_ERR;
   for (i = 1; i <= n; i++) {
      if ((err = mp_init(&(s_schoenhagecache[i]))) != MP_OKAY) {
         return err;
      }
      if ((err = mp_sqr(&(s_schoenhagecache[i - 1]),  &(s_schoenhagecache[i]))) != MP_OKAY) {
         return err;
      }
   }
   s_schoenhagecache_len = i;


   /* we need a defined starting point */
   string = (char *)s_mp_memset(string, '\0', maxlen);
   sign = a->sign;
   if (sign == MP_NEG) {
      *string = '-';
      string++;
      *string = '\0';
   }
   a_bis = *a;
   a_bis.sign = MP_ZPOS;

   if ((err = s_mp_get_str_intern(&a_bis, string, 0, base, buffer_size,
                                  s_schoenhagecache)) != MP_OKAY)   goto LBL_ERR;

   if (written != NULL) {
      *written = (size_t)s_mp_strlen(string) + 1;
   }
LBL_ERR:
   for (i = 0; i < s_schoenhagecache_len; i++) {
      mp_clear(&(s_schoenhagecache[i]));
   }
   return err;
}

#endif
