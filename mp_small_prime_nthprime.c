#include "tommath_private.h"
#ifdef MP_SMALL_PRIME_NTHPRIME_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

static mp_digit s_mp_bitlength(mp_digit n)
{
   mp_digit bits, x = n;
   for (bits = 0; x != 0; bits++) {
      x >>= 1;
   }
   return bits;
}

mp_err mp_small_prime_nthprime(mp_digit n, mp_digit *d)
{
   mp_err err = MP_OKAY;
   mp_digit upper_bound, low, mid, high, res, ilog2, iloglog2, primecount;

   if (n >= PRIMECOUNT_MAX) {
      return MP_VAL;
   }

   if (n < 10) {
      upper_bound = 50;
   } else {
      /* Dusart's bound: n * (log(n) + log(log(n))) -1 for all n > 6 */
      ilog2 = s_mp_bitlength(n);
      iloglog2 = s_mp_bitlength(ilog2);

      /* 1143/1649 ~ 0.69314736 and log(2) ~ 0.69314718, 6 decimals are close enough
         49180508/70952475 for 16 decimals.
         There is also 710/2^10 for three decimals, 726818/2^20 for six decimals, and
         1488522237/2^32 for ten decimals.
      */
      /* Does not look pretty but works and is fast enough, bottleneck is the primecounter() */
      upper_bound = n * (((ilog2 / 1649) * 1143) + ((iloglog2 / 1649) * 1143)) ;
      /* Underflow, try again */
      if (upper_bound < n) {
         upper_bound = n * (((ilog2 * 1143) / 1649) + ((iloglog2 * 1143) / 1649)) ;
      }
      /* overflow, n too large, just set to max */
      if (upper_bound < n) {
         upper_bound =  PRIMECOUNT_MAX;
      }
   }
   low = 2u;
   mid = 0;
   high = upper_bound;
   res = high;

   while (low <= high) {
      mid = low + (high - low) / 2;
      if ((err = mp_small_prime_primecount(mid, &primecount)) != MP_OKAY) {
         return err;
      }
      if (primecount >= n) {
         res = mid;
         high = mid - 1u;
      } else {
         low = mid + 1;
      }
   }
   *d = res;
   return MP_OKAY;
}

#endif
