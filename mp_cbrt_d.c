#include "tommath_private.h"
#ifdef MP_CBRT_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

void mp_cbrt_d(const mp_digit n, mp_digit *r)
{
   int bits = 0;
   mp_digit x, x_sq, x_next;

   if (n == 0u) {
      *r = 0u;
      return;
   }

   if (n < 8u) {
      *r = 1u;
      return;
   }

   bits = s_mp_count_bits_d(n);
   x = (mp_word)1u << ((bits + 2u) / 3u);

   while (true) {
      x_sq = x * x;
      x_next = (2u * x + n / x_sq) / 3u;
      if (x_next >= x) {
         if (x * x * x > n) {
            *r = x - 1u;
            return;
         }
         *r = x;
         return;
      }
      x = x_next;
   }
   /* Won't reach */
   *r = MP_MASK;
   return;
}



#endif
