#include "tommath_private.h"
#ifdef S_MP_MULMOD_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



mp_digit s_mp_mulmod_d(mp_digit a, mp_digit b, mp_digit m)
{
   return (mp_digit)(((mp_word)a * (mp_word)b) % (mp_word)m);
}



#endif
