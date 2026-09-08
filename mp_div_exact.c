#include "tommath_private.h"
#ifdef MP_DIV_EXACT_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */


/*
   Compute the modular inverse b^(-1) mod 2^MP_DIGIT_BIT
   for a mp_digit with Hensel lifting (because 2 is a prime number).
*/
static mp_digit s_mp_invmod_d(mp_digit b)
{
   mp_digit x = b;
   x = (x * (2 - b * x)) & MP_MASK;
   x = (x * (2 - b * x)) & MP_MASK;
   x = (x * (2 - b * x)) & MP_MASK;
   /* MP_xBIT with x in 31,32 */
#if MP_DIGIT_BIT > 16
   x = (x * (2 - b * x)) & MP_MASK;
#endif
   /* MP_64BIT */
#if MP_DIGIT_BIT > 32
   x = (x * (2 - b * x)) & MP_MASK;
#endif
   return x;
}

/*
    Exact integer division
    Computes Q = A / B using Jebelean's algorithm

    @inproceedings{krandick1994bidirectional,
       title={Bidirectional Exact Integer Division.},
       author={Krandick, Werner and Jebelean, Tudor},
       booktitle={PASCO},
       pages={264--272},
       year={1994},
       organization={World Scientific}
    }

    @article{jebelean1993algorithm,
       title={An algorithm for exact division},
       author={Jebelean, Tudor},
       journal={Journal of symbolic computation},
       volume={15},
       number={2},
       pages={169--180},
       year={1993},
       publisher={Elsevier}
    }
 */

mp_err mp_div_exact(const mp_int *A, const mp_int *B, mp_int *Q)
{
   mp_err err;
   mp_int A_scaled, B_scaled, tmp;
   int i, j, k, trailing_zeros = 0, size_q = 0;
   mp_digit b0_inv, q_limb, prod_limb, carry_limb;
   mp_word carry, prod, maskplusone;

   if (mp_iszero(B)) {
      return MP_VAL;
   }

   if ((err = mp_init_multi(&A_scaled, &B_scaled, &tmp, NULL)) != MP_OKAY) {
      return err;
   }

   trailing_zeros = mp_cnt_lsb(B);
   if (trailing_zeros > 0) {
      /* Division is exact, so "A" has, at least, the same amount of trailing_zeros */
      if ((err = mp_div_2d(A, trailing_zeros, &A_scaled, NULL)) != MP_OKAY)                             goto LTM_ERR;
      if ((err = mp_div_2d(B, trailing_zeros, &B_scaled, NULL)) != MP_OKAY)                             goto LTM_ERR;
   } else {
      if ((err = mp_copy(A, &A_scaled)) != MP_OKAY)                                                     goto LTM_ERR;
      if ((err = mp_copy(B, &B_scaled)) != MP_OKAY)                                                     goto LTM_ERR;
   }

   size_q = A_scaled.used - B_scaled.used + 1;
   if (size_q <= 0) {
      mp_zero(Q);
      err = MP_VAL;
      goto LTM_ERR;
   }
   if ((err = mp_grow(Q, size_q)) != MP_OKAY)                                                            goto LTM_ERR;
   Q->used = size_q;

   if ((err = mp_copy(&A_scaled, &tmp)) != MP_OKAY)                                                      goto LTM_ERR;

   b0_inv = s_mp_invmod_d(B_scaled.dp[0]);
   /* What is cheaper, an addition or a shift? */
   maskplusone = (mp_word)1 << MP_DIGIT_BIT;
   for (i = 0; i < size_q; i++) {
      /* Current quotient limb */
      q_limb = (tmp.dp[i] * b0_inv) & MP_MASK;
      Q->dp[i] = q_limb;

      /* tmp = tmp - q_limb * B_scaled */
      carry = 0;
      for (j = 0; j < B_scaled.used; j++) {
         if (i + j >= tmp.used) {
            break;
         }
         /* Handle borrows carfully to avoid trouble with negative results
            which would be a catastrophic failure with unsigned integers here */
         prod = (mp_word)q_limb * (mp_word)B_scaled.dp[j] + carry;
         prod_limb = (mp_digit)(prod & MP_MASK);
         carry = prod >> MP_DIGIT_BIT;
         if (tmp.dp[i + j] < prod_limb) {
            tmp.dp[i + j] = (mp_digit)((tmp.dp[i + j] + maskplusone - prod_limb) & MP_MASK);
            carry++;
         } else {
            tmp.dp[i + j] = (mp_digit)((tmp.dp[i + j] - prod_limb) & MP_MASK);
         }
      }

      /* Propagate those carries/borrows. Carefully. */
      k = i + B_scaled.used;
      while (carry > 0 && k < tmp.used) {
         carry_limb = (mp_digit)(carry & MP_MASK);
         if (tmp.dp[k] < carry_limb) {
            tmp.dp[k] = (mp_digit)((tmp.dp[k] + maskplusone - carry_limb) & MP_MASK);
            carry = (carry >> MP_DIGIT_BIT) + 1;
         } else {
            tmp.dp[k] = (mp_digit)((tmp.dp[k] - carry_limb) & MP_MASK);
            carry = carry >> MP_DIGIT_BIT;
         }
         k++;
      }
   }

   for (k = size_q; k < tmp.used; k++) {
      /* Remainder != 0 */
      if (tmp.dp[k] != 0) {
         /* TODO: add a MP_INEXACT error? */
         err = MP_VAL;
         goto LTM_ERR;
      }
   }

   mp_clamp(Q);

LTM_ERR:
   mp_clear_multi(&A_scaled, &B_scaled, &tmp, NULL);
   return err;
}













#endif
