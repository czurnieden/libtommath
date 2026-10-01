#include "tommath_private.h"
#ifdef S_MP_POWMOD_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

mp_digit s_mp_powmod_d(mp_digit b, mp_digit e, mp_digit modulus)
{
   mp_digit res = 1;

   b %= modulus;
   while (e > 0) {
      res %= modulus;
      if (e & 1) {
         res = s_mp_mulmod_d(res, b, modulus);
      }
      b %= modulus;
      b = s_mp_mulmod_d(b, b, modulus);
      e >>= 1;
   }
   res %= modulus;
   return res;
}

#endif
