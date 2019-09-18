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

/* There is a possibility of a bit of tunability here, I may note. */
#define SCHOENHAGE_CONVERSION_CUT  10

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
                                  char *ls, mp_int *s_schoenhagecache)
{
   int b, n, i;
   int ed;
   mp_err err;
   mp_int q, r;
   char *_s;

   /* Use naive method for the leaves */
   if (a->used <= SCHOENHAGE_CONVERSION_CUT) {
      ls = s_mp_memset(ls,'\0',SCHOENHAGE_CONVERSION_CUT * MP_DIGIT_BIT + 1);
      /*
          TODO: replace mp_to_radix with an internal version when the new
                version comes out that has all buffer allocation done
                internally.

                Or refactor (just get rid of the "ls" char-array) to use
                the new one. That would also remove the need for s_mp_memset
                and one s_mp_strlen call.
       */
      if ((err = mp_to_radix(a, ls, SIZE_MAX, base)) != MP_OKAY) {
         return err;
      }
      /*
          I think that there is a slight chance to put those leading zeros
          in their places in a more elegant way ;-)
       */
      i = s_mp_strlen(ls);
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
      string = s_mp_strncat(string, ls, SIZE_MAX);
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
                                  base, ls, s_schoenhagecache)) != MP_OKAY)    goto LBL_ERR;
   if ((err = s_mp_get_str_intern(&r, string, ed,
                                  base, ls, s_schoenhagecache)) != MP_OKAY)             goto LBL_ERR;

LBL_ERR:
   mp_clear_multi(&q, &r, NULL);
   return MP_OKAY;
}

/* The buffer "string" must be allocated in advance */
/* TODO: do it internally when the new mp_to_radix comes out */
mp_err mp_get_str(const mp_int *a, char *string, int base)
{
   mp_sign sign;
   mp_err e = MP_OKAY;
   mp_int a_bis;
   /* maximum from base 2 */
   /* mp_radix_size was too slow for it at that time */
   char s[SCHOENHAGE_CONVERSION_CUT * MP_DIGIT_BIT + 1];

   mp_int s_schoenhagecache[(sizeof(int) * CHAR_BIT) * sizeof(mp_int)];
   int s_schoenhagecache_len = 0, i, n, b;

   /* check range of the maxlen, radix */
   if ((base < 2) || (base > 64)) {
      return MP_VAL;
   }

   /* quick out if its zero */
   if (MP_IS_ZERO(a)) {
      *string++ = '0';
      *string = '\0';
      return MP_OKAY;
   }

   b = mp_count_bits(a);
   if (b <= (SCHOENHAGE_CONVERSION_CUT * MP_DIGIT_BIT)) {
      return  mp_to_radix(a, string, SIZE_MAX, base);
   }

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

       TODO: we don't need the entries of size < (SCHOENHAGE_CONVERSION_CUT * MP_DIGIT_BIT)
             (in their respective bases, of course)

       TODO: The entries get used multiple times, so it is a good idea to replace the entries
             with their approximate reciprocals and use Barrett division.

    */
   if ((e = mp_init_set(&(s_schoenhagecache[0]), (mp_digit)(base))) != MP_OKAY) goto LBL_ERR;
   for (i = 1; i <= n; i++) {
      if ((e = mp_init(&(s_schoenhagecache[i]))) != MP_OKAY) {
         return e;
      }
      if ((e = mp_sqr(&(s_schoenhagecache[i - 1]),  &(s_schoenhagecache[i]))) != MP_OKAY) {
         return e;
      }
   }
   s_schoenhagecache_len = i;


   /* we need a defined starting point */
   *string = '\0';
   sign = a->sign;
   if (sign == MP_NEG) {
      *string = '-';
      string++;
      *string = '\0';
   }
   a_bis = *a;
   a_bis.sign = MP_ZPOS;

   if ((e = s_mp_get_str_intern(&a_bis, string, 0, base, s,
                                s_schoenhagecache)) != MP_OKAY)   goto LBL_ERR;


LBL_ERR:
   for (i = 0; i < s_schoenhagecache_len; i++) {
      mp_clear(&(s_schoenhagecache[i]));
   }
   return e;
}

#endif
