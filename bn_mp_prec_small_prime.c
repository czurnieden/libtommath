#include "tommath_private.h"
#ifdef BN_MP_PREC_SMALL_PRIME_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

/*
 * Mimics behaviour of Pari/GP's precprime(n)
 * If n is prime set *result to n else set *result to first prime < n
 * and 0 in case of error
 */
mp_err mp_prec_small_prime(mp_sieve_prime n, mp_sieve_prime *result, mp_sieve *sieve)
{
   mp_sieve_prime ret = 0;
   int e = MP_OKAY;

   if (n == 2) {
      *result = 2;
      return e;
   }

   if (n < 2) {
      *result = 0;
      return e;
   }

   for (; ret == 0; n--) {
      if ((e = mp_is_small_prime(n, &ret, sieve)) != MP_OKAY) {
         *result = 0;
         return e;
      }
   }
   *result = n + 1;

   return e;
}
#endif

