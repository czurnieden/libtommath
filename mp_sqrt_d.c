#include "tommath_private.h"
#ifdef MP_SQRT_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

/* Just measure, fast enough */
#define MP_DIGIT_BITS (sizeof(mp_digit)*CHAR_BIT)
mp_err mp_sqrt_d(const mp_digit n, mp_digit *r)
{
   mp_err err = MP_OKAY;
   mp_digit s, rem, root;

   if (r == NULL) {
      return MP_VAL;
   }

   if (n <= 1u) {
      *r = n;
      return err;
   }

   /* highest power of four <= n */
   s = (mp_digit)(1ull << (MP_DIGIT_BITS - 2));
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
