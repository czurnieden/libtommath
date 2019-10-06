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

/* Tune-able */
/* Size in bits */
#define BARRETT_CUT 1000


/* based on SCHOENHAGE_CONVERSION_CUT=600 */
/*
   Computed with (Pari/GP)
      len(cut, radix) = ceil(log(2^(cut+1))/log(radix)) + 1
      for(i = 2, 64, k = len(600,i);printf(k ", "))

    Although the table is hardcoded for 600 bit a change to SCHOENHAGE_CONVERSION_CUT
    is directly proportional to the sizes in the table, e.g.: a raise of
    the SCHOENHAGE_CONVERSION_CUT value by 10% gives rise to the values in the table by
    10%, too.
    Adding some angst-allowance for the unavoidable rounding errors is highly recommended.
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

/* Cache */
typedef struct s_mp_schoenhage {
   mp_int entry;
   mp_int reciprocal;
   int shift_value;
   mp_bool is_reciprocal;
} s_mp_schoenhage;

static mp_err s_mp_cache_entry_init(s_mp_schoenhage *entry)
{
   mp_err err = MP_OKAY;

   mp_init(&(entry->entry));
   mp_init(&(entry->reciprocal));
   entry->shift_value = 0;
   entry->is_reciprocal = MP_NO;
   return err;
}

static void s_mp_cache_entry_clear(s_mp_schoenhage *entry)
{
   mp_clear(&(entry->entry));
   mp_clear(&(entry->reciprocal));
}


/* Barrett division as described in Brent/Zimmernman "Modern Computer Arithmetic" p. 64f */
/* Computes reciprocal c = floor((2^(k))/a) for Barrett division */
static mp_err s_mp_barrett_reciprocal(const mp_int *a, int beta2, mp_int *c, int *shift_value)
{
   mp_err err = MP_OKAY;
   mp_int t, c1;
   int shift;

   mp_init_multi(&t, &c1, NULL);
   /* I = floor(beta^2/B) */

   mp_2expt(&t, beta2);
   /* 
     This is a very naive approach but is already a magnitude faster than mp_to_radix.
     You'll get another speed-up with fast division, obviously, but only if the denominator
     is larger than the Karatsuba cut-off.

     But there are methods to speed  that 2^k/b up a bit e.g.: with Newton's method to
     get an reciprocal.

     Example:

        static mp_err s_mp_newton_inverse(const mp_int *a, int exp, mp_int *inverse, int *shift_value){
           mp_err err = MP_OKAY;
           mp_int xn, twokp1, t1, t2;
           int start_value, i = 0;
        
           
           mp_init_multi(&xn, &t1, &t2, &twokp1, NULL);
        
           start_value = exp - mp_count_bits(a);
           
           mp_2expt(&xn, start_value);
           mp_2expt(&twokp1, exp + 1);
        
           for (;;) {
              mp_copy(&xn, &t1);
              mp_mul(&xn, a, &t2);
              mp_sub(&twokp1, &t2, &t2);
              mp_mul(&xn, &t2, &xn);
              mp_div_2d(&xn, exp, &xn, NULL);
              mp_sub(&xn, &t1, &t2);
              t2.sign = MP_ZPOS;
              if(mp_cmp_d(&t2, 2) == MP_LT) {
                 break;
              }
           }
          mp_exch(&xn, inverse);
          *shift_value = exp;
        
          mp_clear_multi(&xn, &t1, &t2, &twokp1,  NULL);
          return err;
        }

     Which can be optimized further e.g.: you don't need the full precision for every step.

     Oh, and yes, we allow for being one-off to avoid getting stuck (Newton's method for integers
     has the habit to get stuck between floor(result) and ceil(result) and loops).
   */
   mp_div(&t, a, &c1, NULL);
 
   *shift_value = beta2;

   mp_exch(&c1, c);
   mp_clear_multi(&t, &c1, NULL);
   return err;
}

/* Barrett division with correction: compute A*I/beta and correct the result */
/* 0 >= A < beta^2  and beta/2 < B < beta */
static mp_err s_mp_barrett_division(const mp_int *A, const mp_int *I, const mp_int *B,
                                    int beta, mp_int *q, mp_int *r)
{
   mp_err err = MP_OKAY;
   /* TODO: use q,r directly */
   mp_int Q, R, A1;

   mp_init_multi(&Q, &R, &A1, NULL);

   /* Q = floor( (A_1 * I) / beta) where A = A_1 * beta + A_0 with 0 <= A_0 < beta */
   mp_mul(A, I, &Q);
   /* The beta here is the same beta as the beta in beta^2 in s_mp_barrett_reciprocal */
   mp_div_2d(&Q, 2 * beta, &Q, NULL);
   /* R = A - Q*B */
   mp_mul(&Q, B, &R);
   mp_sub(A, &R, &R);
   /*
      If done correctly, this loop does not run more than three times.
      The input is not always *exactly* in the necessary ranges but the number
      of iterations should be in the single digits.
    */
   /* while R >= B */
   while (mp_cmp(&R, B) != MP_LT) {
      /* Q = Q + 1 */
      mp_incr(&Q);
      /* R = R - B */
      mp_sub(&R, B, &R);
   }
   mp_exch(&Q, q);
   mp_exch(&R, r);
   mp_clear_multi(&Q, &R, &A1, NULL);
   return err;
}

static mp_err s_mp_get_str_intern(mp_int *a, char *string, int digits, int base,
                                  size_t size, s_mp_schoenhage *s_schoenhagecache,
                                  mp_bool use_barrett)
{
   int b, beta, n, i, shift;
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
   n = s_ilog2(b / (2 * log_table[base]));
   if (n >= (int)(sizeof(int) * CHAR_BIT) - 1) {
      return MP_VAL;
   }

   if ((err = mp_init_multi(&q, &r, NULL)) != MP_OKAY) {
      return err;
   }
   /*
       Divide a chunk of the input by the correct basepower from the cache.
       This is the point where a fast division algorithm does come handy.
    */
   if ((use_barrett == MP_YES)) {
      if (s_schoenhagecache[n].is_reciprocal == MP_NO) {
         b = (b * 3) / 2;
         beta = (b)/ 2;
         /* We need b = 2*beta an 2^b > a*/
         if ((2*beta) < b) {
            beta++;
         }
         s_mp_barrett_reciprocal(&(s_schoenhagecache[n].entry),
                                 2*beta,
                                 &(s_schoenhagecache[n].reciprocal), &shift);
         s_schoenhagecache[n].shift_value = shift/2;
         s_schoenhagecache[n].is_reciprocal = MP_YES;
      }
      s_mp_barrett_division(a, &(s_schoenhagecache[n].reciprocal),
                            &(s_schoenhagecache[n].entry),
                            s_schoenhagecache[n].shift_value, &q, &r);
   } else {
      if ((err = mp_div(a, &(s_schoenhagecache[n].entry), &q, &r)) != MP_OKAY)       goto LBL_ERR;
   }
   ed = 1 << n;
   /* Rinse and repeat with quotient and remainder */
   if ((err = s_mp_get_str_intern(&q, string, digits - ed,
                                  base, size, s_schoenhagecache,
                                  use_barrett)) != MP_OKAY)                 goto LBL_ERR;
   if ((err = s_mp_get_str_intern(&r, string, ed,
                                  base, size, s_schoenhagecache,
                                  use_barrett)) != MP_OKAY)                 goto LBL_ERR;

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
   mp_bool use_barrett = MP_NO;
   mp_int a_bis;

   int s_schoenhagecache_len = 0, i, n, b;
   size_t buffer_size, wrote = 0u;

   s_mp_schoenhage s_schoenhagecache[(sizeof(int) * CHAR_BIT) * sizeof(s_mp_schoenhage)];

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

   /* Compute size of buffer */
   buffer_size = (buf_size[base] * SCHOENHAGE_CONVERSION_CUT) / SCHOENHAGE_CONVERSION_CUT_BASE;

   /* Compute the max exponent for the base necessary */
   n = s_ilog2(b / (2 * log_table[base]));
   if (n >= (int)(sizeof(int) * CHAR_BIT) - 1) {
      return MP_VAL;
   }
   /*
       Fill the cache with squares of squares of "base". Example with base 10 (ten)

       s_schoenhagecache[0] = 10
       s_schoenhagecache[1] = s_schoenhagecache[0]^2 = 100
       s_schoenhagecache[2] = s_schoenhagecache[1]^2 = 10000
       ...

       TODO: we don't need the smaller entries (in their respective bases, of course),
             they are not used.

       A Q&D benchmark gave some values around half a million to a million bits as a cut-off
       to use Barrett-divisions but YMMV, as always, so it needs a tuning mechanism.
    */
   use_barrett = (b >= BARRETT_CUT)? MP_YES: MP_NO;

   s_mp_cache_entry_init(&(s_schoenhagecache[0]));
   mp_set(&(s_schoenhagecache[0].entry), (mp_digit)(base));
   s_schoenhagecache[0].is_reciprocal = MP_NO;

   for (i = 1; i <= n; i++) {
      if ((err = s_mp_cache_entry_init(&(s_schoenhagecache[i]))) != MP_OKAY) {
         return err;
      }
      if ((err = mp_sqr(&(s_schoenhagecache[i - 1].entry),  &(s_schoenhagecache[i].entry))) != MP_OKAY) {
         return err;
      }
      s_schoenhagecache[n].is_reciprocal = MP_NO;
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
                                  s_schoenhagecache, use_barrett)) != MP_OKAY)   goto LBL_ERR;
   if (written != NULL) {
      *written = (size_t)s_mp_strlen(string) + 1;
   }
LBL_ERR:
   for (i = 0; i < s_schoenhagecache_len; i++) {
      s_mp_cache_entry_clear(&(s_schoenhagecache[i]));
   }
   return err;
}

#endif
