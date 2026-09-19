#include "tommath_private.h"
#ifdef S_MP_COUNT_BITS_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



int s_mp_count_bits_d(mp_digit n)
{
   int bits = 0;
   while (n > 0) {
      bits++;
      n >>= 1;
   }
   return bits;
}


#endif
