#include "tommath_private.h"
#ifdef MP_SMALL_PRIME_PRIMECOUNT_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

mp_err mp_small_prime_primecount(mp_digit n, mp_digit *d)
{
   mp_err err = MP_OKAY;
   mp_digit sqrt_n, m = 0u, v, p, p2, sp_prev, vp, idx_vp, i = 0u, j;
   mp_digit *V, *S;

   if (n < 2) {
      *d = 0u;
      return err;
   }

   /* Really necessary? */
   if (n > (MP_MASK - 3)) {
      n = n - 3;
   }

   if ((err = mp_sqrt_d(n, &sqrt_n)) != MP_OKAY) {
      return err;
   }
   V = (mp_digit *)MP_MALLOC(sizeof(mp_digit) * (2u * sqrt_n + 2u));
   if (V == NULL) {
      return MP_MEM;
   }

   for (i = 1u; i <= n; i = j + 1u) {
      v = n / i;
      V[m++] = v;
      j = n / v;
   }

   S = (mp_digit *)MP_MALLOC(sizeof(mp_digit) * m);
   if (S == NULL) {
      MP_FREE(V, sizeof(mp_digit) * (2u * sqrt_n + 2u));
      return MP_MEM;
   }

   for (i = 0u; i < m; i++) {
      S[i] = V[i] - 1;
   }

   for (p = 2u; p <= sqrt_n; p++) {
      if (S[m - p] > S[m - (p - 1u)]) {
         sp_prev = S[m - (p - 1u)];
         p2 = p * p;
         for (i = 0u; i < m; i++) {
            if (V[i] < p2) {
               break;
            }
            vp = V[i] / p;
            if (vp <= sqrt_n) {
               idx_vp = m - vp;
            } else {
               idx_vp = n / vp - 1u;
            }
            if (vp > sqrt_n) {
               idx_vp = n / vp - 1;
            } else {
               idx_vp = m - vp;
            }
            S[i] -= (S[idx_vp] - sp_prev);
         }
      }
   }
   *d = S[0];
   MP_FREE(V, sizeof(mp_digit) * (2u * sqrt_n + 2u));
   MP_FREE(S, sizeof(mp_digit) * m);
   return err;
}

#endif
