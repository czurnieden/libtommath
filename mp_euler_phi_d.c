#include "tommath_private.h"
#ifdef MP_EULER_PHI_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */




/* Compute Euler's totient function.
   Optimized for random input using only little memory.

   For sequential uses over a larger range other algorithms
   are better suited (A Sundaram sieve (sort of, kinda) together
   with Euler's product formula). Although they use a lot of memory
   for larger ranges they are quite fast.

   Actual factoring as done here is quite slow.
 */

/* TODO: use own sieve or one given? */
mp_err mp_euler_phi_d(const mp_digit x, mp_digit *res, mp_erat_sieve *sieve)
{
   mp_digit result = x, x2 = x;
   ERAT_UINT p = 2u;
   mp_digit factor;
   mp_err err = MP_OKAY;

   if (x <= 1) {
      *res = x;
      return err;
   }

   if (s_mp_isprime_d(x)) {
      *res = x- 1u;
      return err;
   }

   /* Use the base sieve only */
   while (p < ERAT_UINT_MAX_SQRT) {
      if (x2 % p == 0) {
         result -= result / p;
         while (x2 % p == 0) {
            x2 /= p;
         }
      }
      if ((err = mp_small_prime_sieve_next_prime(p + 1u, &p, sieve)) != MP_OKAY) {
         if (err != MP_OVF) {
            break;
         }
         goto LTM_ERR;
      }
   }

   if (x2 == 1) {
      *res = result;
      return err;
   }

   if (s_mp_isprime_d(x2)) {
      result -= result / x2;
      goto LTM_END;
   }

   while (true) {
      /* get larger factor*/
      factor = s_mp_pollard_rho_d((mp_digit)x2, 1000000);
      if (factor == 0) {
         goto LTM_ERR;
      }
      /* One time only, the small primes are sieved out already */
      if (!s_mp_isprime_d(factor)) {
         factor = s_mp_pollard_rho_d((mp_digit)x2, 1000000);
         if (factor == 0) {
            goto LTM_ERR;
         }
      }

      result -= result / factor;
      while (x2 % factor == 0) {
         x2 /= factor;
      }
      if ((x2 == 1) || (s_mp_isprime_d(x2))) {
         break;
      }
   }
   if (s_mp_isprime_d(x2)) {
      result -= result / x2;
   }

LTM_END:
   *res = result;
   return err;
LTM_ERR:
   return err;
}



#endif
