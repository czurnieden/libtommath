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

/* There is a possibility of a bit of tunability here, I may note. */
#define SCHOENHAGE_CONVERSION_CUT  10

/*
    Threadsafety was not the most urgent thing for this prototype.
    You need one, though, the difference is significant for even moderately
    large numbers but doing it in the root (mp_get_str) would suffice.
 */
static mp_int *s_schoenhagecache;
static int s_schoenhagecache_len;
static int s_schoenhagecache_base;

static void s_free_schoenhage_cache(void)
{
   int i;
   for (i = 0; i < s_schoenhagecache_len; i++)
      mp_clear(&(s_schoenhagecache[i]));
   free(s_schoenhagecache);
}

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

static mp_err s_mp_get_str_intern(mp_int *a, char *string, int digits, int base, char *ls)
{
   int b, n, i;
   int ed;
   mp_err err;
   mp_int q, r;

   /* Use naive method for the leaves */
   if (a->used <= SCHOENHAGE_CONVERSION_CUT) {
      ls = memset(ls,'\0',SCHOENHAGE_CONVERSION_CUT * MP_DIGIT_BIT + 1);
      if ((err = mp_to_radix(a, ls, SIZE_MAX, base)) != MP_OKAY) {
         return err;
      }
      /*
          I think that there is a slight chance to put those leading zeros
          in their places in a more elegant way ;-)
       */
      if ((strlen(ls) < (unsigned) digits) && (strlen(string) > 0)) {
         for (i = strlen(ls); i < digits; i++) {
            string = strncat(string, "0", 1);
         }
      }
      string = strcat(string, ls);
      return MP_OKAY;
   }

   b = mp_count_bits(a);
   /* Compute the exponent for the base necessary */
   n = s_ilog2(b / (2 * log_table[base])) - 1;
   if (n >= (int)(sizeof(int) * CHAR_BIT) - 1) {
      return MP_VAL;
   }

   /* Build the cache if it is empty. */
   if (s_schoenhagecache_len == 0) {
      /*
          There is no need for a malloc, the size is known at compile time
          although it is not very small (e.g.: 768 bytes for a 64 bit architecture
          with 4-byte "int"'s, 8-byte pointers and a 4-byte padding in "mp_int")
          so malloc might be useful where memory is scarce.
       */
      s_schoenhagecache = malloc((sizeof(int) * CHAR_BIT) * sizeof(mp_int));
      if (s_schoenhagecache == NULL) {
         return MP_MEM;
      }
      /*
          Fill the cache with squares of squares of "base". Example with base 10 (ten)

          s_schoenhagecache[0] = 10
          s_schoenhagecache[1] = s_schoenhagecache[0]^2 = 100
          s_schoenhagecache[2] = s_schoenhagecache[1]^2 = 10000
          ...
       */
      if( (err = mp_init_set(&(s_schoenhagecache[0]), (mp_digit)(base)) ) != MP_OKAY) goto LBL_ERR;
      for (i = 1; i < n; i++) {
         if ((err = mp_init(&(s_schoenhagecache[i]))) != MP_OKAY) {
            return err;
         }
         if ((err =
                 mp_sqr(&(s_schoenhagecache[i - 1]),
                        &(s_schoenhagecache[i]))) != MP_OKAY) {
            return err;
         }
      }
      s_schoenhagecache_len = i;
   }
   /* Cache might be insuffiently filled (e.g.: from last time converting a smaller number) */
   if (n >= s_schoenhagecache_len) {
      for (i = s_schoenhagecache_len; i <= n; i++) {
         if ((err = mp_init(&(s_schoenhagecache[i]))) != MP_OKAY) {
            return err;
         }
         if ((err =
                 mp_sqr(&(s_schoenhagecache[i - 1]),
                        &(s_schoenhagecache[i]))) != MP_OKAY) {
            return err;
         }
      }
      s_schoenhagecache_len = i;
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
   if ((err = s_mp_get_str_intern(&q, string, digits - ed, base, ls)) != MP_OKAY)    goto LBL_ERR;
   if ((err = s_mp_get_str_intern(&r, string, ed, base, ls)) != MP_OKAY)             goto LBL_ERR;

LBL_ERR:
   mp_clear_multi(&q, &r, NULL);
   return MP_OKAY;
}

/* The buffer "string" must be sufficiently allocated in advance */
mp_err mp_get_str(const mp_int *a, char *string, int base)
{
   mp_sign sign;
   mp_err e;
   mp_int a_bis;
   /* maximum from base 2 */
   /* mp_radix_size was too slow for it at that time */
   char s[SCHOENHAGE_CONVERSION_CUT * MP_DIGIT_BIT + 1];

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

   /* made static for simplicity */
   /*if (s_schoenhagecache_base != base) {*/
      s_free_schoenhage_cache();
   /*}*/
   s_schoenhagecache_base = base;

   if ((e = s_mp_get_str_intern(&a_bis, string, 0, base, s)) != MP_OKAY) {
      return e;
   }

   return MP_OKAY;
}




#endif
