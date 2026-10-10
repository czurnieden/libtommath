#include "tommath_private.h"
#ifdef S_MP_CNT_LSB_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



int s_mp_cnt_lsb_d(mp_digit v)
{
   int count = 0;
   if (v == 0) {
      return (sizeof(mp_digit) * CHAR_BIT);
   }
   while ((v & 1u) == 0) {
      count++;
      v >>= 1;
   }
   return count;
}



#endif
