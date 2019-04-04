#include "tommath_private.h"
#ifdef BN_MP_NEXT_SMALL_PRIME_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

/*
 * Mimics behaviour of Pari/GP's nextprime(n)
 * If n is prime set *result to n else set *result to first prime > n
 * and 0 in case of error
 */
mp_err mp_next_small_prime(mp_sieve_prime n, mp_sieve_prime *result, mp_sieve *sieve)
{
   mp_sieve_prime ret = 0;
   int e = MP_OKAY;

   if (n < 2) {
      *result = 2;
      return e;
   }

   for (; ret == 0 && n != 0; n++) {
      if (n > MP_SIEVE_BIGGEST_PRIME) {
         return MP_SIEVE_MAX_REACHED;
      }
      /* just call mp_is_small_prime(), it does all of the heavy work */
      if ((e = mp_is_small_prime(n, &ret, sieve)) != MP_OKAY) {
         *result = 0;
         return e;
      }
   }
   *result = n - 1;
   return e;
}
#endif

