#include "tommath_private.h"
#ifdef S_MP_GCD_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



mp_digit s_mp_gcd_d(mp_digit a, mp_digit b)
{
   mp_digit x = a;
   mp_digit y = b;
   mp_digit g = 1u, t;

   if (a == b) {
      return a;
   }

   if ((a == 0u) && (b != 0u)) {
      return b;
   }

   if ((a != 0u) && (b == 0u)) {
      return a;
   }
   if ((a == 0u) && (b == 0u)) {
      return 0;
   }

   while (!(x & 1u) && !(y & 1u)) {
      x = x >> 1u;
      y = y >> 1u;
      g = g << 1u;
   }
   while (x != 0u) {
      while (!(x & 1u)) {
         x >>= 1u;
      }
      while (!(y & 1u)) {
         y >>= 1u;
      }
      t = (y > x) ? (y - x) >> 1u : (x - y) >> 1u ;
      if (x >= y) {
         x = t;
      } else y = t;
   }
   return (g * y);
}


#endif
