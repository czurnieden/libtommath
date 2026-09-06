#include "tommath_private.h"
#ifdef MP_VALUATION_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

mp_err mp_valuation_d(const mp_int *a, const mp_digit p, mp_digit *r)
{
   mp_err err = MP_OKAY;
   mp_int tmp;
   mp_digit k = 0u, rem = 0u;

   if ((p < 2) || mp_iszero(a)) {
      *r = mp_cnt_lsb(a);
      return err;
   }
   if ((err = mp_init_copy(&tmp, a)) != MP_OKAY) {
      return err;
   }
   while (true) {
      if ((err = mp_div_d(&tmp, p, &tmp, &rem)) != MP_OKAY)                                               goto LTM_ERR;
      if (rem == 0u) {
         k++;
      } else {
         break;
      }
   }
   *r = k;

LTM_ERR:
   mp_clear(&tmp);
   return err;
}

#endif
