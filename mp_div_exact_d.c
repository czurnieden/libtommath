#include "tommath_private.h"
#ifdef MP_DIV_EXACT_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#define DEBUG_PRINT(ERROR_NUMBER, ERROR_GOTO)\
        do{\
           fprintf(stderr, "%s %d in %s: %s\n",\
           __FILE__, __LINE__, __FUNCTION__,\
            mp_error_to_string((ERROR_NUMBER)));\
           goto ERROR_GOTO;\
        }while(0)


static mp_err s_mp_div_exact_odd_d(const mp_int *a, mp_digit d, mp_int *q)
{
   mp_err err = MP_OKAY;
   mp_digit d_inv, limb_a, limb_q;
   mp_word carry, prod;
   int i;

   if ((err = mp_grow(q, a->used)) != MP_OKAY)                                                           DEBUG_PRINT(err,
            LTM_ERR);
   q->used = a->used;

   d_inv = s_mp_invmod_d(d);
   carry = 0;

   for (i = 0; i < a->used; i++) {
      limb_a = a->dp[i];

      /* (limb_q * d + carry) = limb_a (mod 2^60) <-> limb_q = (limb_a - carry) * d_inv (mod 2^60) */
      if (limb_a >= (carry & MP_MASK)) {
         limb_q = (mp_digit)((((mp_word)limb_a - (carry & MP_MASK)) * d_inv) & MP_MASK);
      } else {
         limb_q = (mp_digit)((((mp_word)limb_a + ((mp_word)1 << MP_DIGIT_BIT) - (carry & MP_MASK)) * d_inv) & MP_MASK);
      }

      q->dp[i] = limb_q;

      prod = ((mp_word)limb_q * (mp_word)d) + carry;

      /* Calculate carry if available */
      if (prod > limb_a) {
         carry = (prod - limb_a) >> MP_DIGIT_BIT;
      } else {
         carry = 0;
      }
   }

   /* Not exact if we have a remainder */
   if (carry != 0) {
      err = MP_VAL;
      DEBUG_PRINT(err, LTM_ERR);
   }
   mp_clamp(q);
LTM_ERR:
   return err;
}

mp_err mp_div_exact_d(const mp_int *a, mp_digit d, mp_int *q)
{
   mp_err err = MP_OKAY;
   int shift_bits = 0;
   mp_digit original_d = d;

   if (d > MP_MASK) {
      err = MP_VAL;
      DEBUG_PRINT(err, LTM_ERR);
   }

   if (d == 0) {
      err = MP_VAL;
      DEBUG_PRINT(err, LTM_ERR);
   }

   if (mp_iszero(a)) {
      mp_zero(q);
      return MP_OKAY;
   }

   /* Shortcut for powers of two. */
   if (MP_IS_2EXPT(d)) {
      shift_bits = s_mp_cnt_lsb_d(d);
      if ((err = mp_div_2d(a, shift_bits, q, NULL)) != MP_OKAY)                                        DEBUG_PRINT(err,
               LTM_ERR);
   }

   if ((d & 1u) == 0) {
      shift_bits = s_mp_cnt_lsb_d(d);
      d >>= shift_bits;
      if ((err = mp_div_2d(a, shift_bits, q, NULL)) != MP_OKAY)                                          DEBUG_PRINT(err,
               LTM_ERR);
      if ((err = s_mp_div_exact_odd_d(q, d, q)) != MP_OKAY)                                              DEBUG_PRINT(err,
               LTM_ERR);
   } else {
      if ((err = s_mp_div_exact_odd_d(a, d, q)) != MP_OKAY)                                              DEBUG_PRINT(err,
               LTM_ERR);
   }

   /* Fallback. Just in case */
   if (err != MP_OKAY) {
      return mp_div_d(a, original_d, q, NULL);
   }
LTM_ERR:
   return err;
}


#endif
