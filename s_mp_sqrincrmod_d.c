#include "tommath_private.h"
#ifdef S_MP_SQRINCRMOD_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



mp_digit s_mp_sqrincrmod_d(mp_digit a, mp_digit m)
{
   return (mp_digit)(((mp_word)a * (mp_word)a + 1u) % (mp_word)m);
}


#endif
