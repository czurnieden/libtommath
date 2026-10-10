#include "tommath_private.h"
#ifdef S_MP_INVMOD_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



/*
   Compute the modular inverse b^(-1) mod 2^MP_DIGIT_BIT
   for a mp_digit with Hensel lifting.
   This works because unsigned overflow is just mod 2^w (w the
   bitsize of the unsigned integer type) which does not change the
   mathematics here.
*/
mp_digit s_mp_invmod_d(mp_digit b)
{
   mp_digit x = b;                   /* 2 */
   x = (x * (2u - b * x)) & MP_MASK; /* 4 */
   x = (x * (2u - b * x)) & MP_MASK; /* 8 */
   x = (x * (2u - b * x)) & MP_MASK; /* 16 */
   /* MP_xBIT with x in 31,32 */
#if MP_DIGIT_BIT > 16
   x = (x * (2u - b * x)) & MP_MASK; /* 32 */
#endif
   /* MP_64BIT */
#if MP_DIGIT_BIT > 32
   x = (x * (2u - b * x)) & MP_MASK; /* 64 */
#endif
   return x;
}



#endif
