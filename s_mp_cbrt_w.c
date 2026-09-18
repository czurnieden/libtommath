#include "tommath_private.h"
#ifdef S_MP_CBRT_W_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



static int mp_count_bits_w(mp_word n)
{
   int bits = 0;
   while (n > 0) {
      bits++;
      n >>= 1;
   }
   return bits;
}

void s_mp_cbrt_w(const mp_word n, mp_word *r)
{
   int bits = 0;
   mp_word x, x_sq, x_next;

   if (n == 0u) {
      *r = 0u;
      return;
   }

   if (n < 8u) {
      *r = 1u;
      return;
   }

   bits = mp_count_bits_w(n);
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
