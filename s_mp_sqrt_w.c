#include "tommath_private.h"
#ifdef S_MP_SQRT_W_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */


/* Just measure, fast enough */
#define MP_DIGIT_BITS (sizeof(mp_word)*CHAR_BIT)
mp_err s_mp_sqrt_w(const mp_word n, mp_word *r)
{
   mp_err err = MP_OKAY;
   mp_word s, rem, root;

   if (r == NULL) {
      return MP_VAL;
   }

   if (n <= 1u) {
      *r = n;
      return err;
   }

   /* highest power of four <= n */
   s = (mp_word)(((mp_word)1u) << (MP_DIGIT_BITS - 2));
   rem = n;
   root = 0u;
   while (s > n) {
      s >>= 2u;
   }
   while (s != 0u) {
      if (rem >= (s | root)) {
         rem -= (s | root);
         root >>= 1u;
         root |= s;
      } else {
         root >>= 1u;
      }
      s >>= 2u;
   }
   *r = root;
   return err;
}

#endif
