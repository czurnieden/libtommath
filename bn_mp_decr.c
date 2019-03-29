#include "tommath_private.h"
#ifdef BN_MP_DECR_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis
 *
 * LibTomMath is a library that provides multiple-precision
 * integer arithmetic as well as number theoretic functionality.
 *
 * The library was designed directly after the MPI library by
 * Michael Fromberger but has been written from scratch with
 * additional optimizations in place.
 *
 * SPDX-License-Identifier: Unlicense
 */

/* Decrement "a" by one like "a--". Changes input! */
#ifdef LTM_USE_EXTRA_FUNCTIONS
int mp_decr(mp_int *a)
{
   int e;
   if (IS_ZERO(a)) {
      mp_set(a,1uL);
      a->sign = MP_NEG;
      return MP_OKAY;
   } else if (a->sign == MP_NEG) {
      a->sign = MP_ZPOS;
      if ((e = mp_incr(a)) != MP_OKAY) {
         a->sign = MP_NEG;
         return e;
      }
      a->sign = MP_NEG;
      return MP_OKAY;
   } else if ((a->used == 1) && (a->dp[0] > 1uL)) {
      a->dp[0]--;
      if (a->dp[0] == 0) {
         mp_zero(a);
      }
      return MP_OKAY;
   } else  {
      return mp_sub_d(a, 1uL,a);
   }
}
#endif
#endif
/* ref:         \$Format:\%D$ */
/* git commit:  \$Format:\%H$ */
/* commit time: \$Format:\%ai$ */
