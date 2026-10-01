#include "tommath_private.h"
#ifdef S_MP_MILLER_RABIN_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

bool s_mp_miller_rabin_d(mp_digit number, mp_digit witness)
{
   mp_digit a = witness, x, y = 1u, d;

   d = number - 1;
   while ((d & 1) == 0) {
      d /= 2;
   }

   x = s_mp_powmod_d(a, d, number);

   if ((x == 1)  || (x == (number - 1))) {
      return true;
   }

   while (d != (number - 1u)) {
      y = s_mp_powmod_d(x, 2u, number);
      d *= 2;

      if ((y == 1u) && (x != 1u) && (x != (number - 1u))) {
         return false;
      }
      x = y;
   }
   if (y != 1u) {
      return false;
   }
   return true;
}

#endif
