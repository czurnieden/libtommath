#include "tommath_private.h"
#ifdef MP_MIN_MULTI_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */


mp_err mp_min_multi(mp_int *out, const mp_int *first, ...)
{
   mp_err err = MP_OKAY;
   mp_int cur_arg;
   mp_int *next_arg;
   va_list args;

   if ((out == NULL) || (first == NULL)) {
      return MP_VAL;
   }
   if ((err = mp_init(&cur_arg)) != MP_OKAY) {
      return err;
   }
   if ((err = mp_copy(first, &cur_arg)) != MP_OKAY)                                                    goto LTM_ERR;
   va_start(args, first);
   next_arg = va_arg(args, mp_int *);
   while (next_arg != NULL) {
      if (mp_cmp(next_arg, &cur_arg) == MP_LT) {
         if ((err = mp_copy(next_arg, &cur_arg)) != MP_OKAY)                                           goto LTM_ERR;
      }
      next_arg = va_arg(args, mp_int *);
   }
   if ((err = mp_copy(&cur_arg, out)) != MP_OKAY)                                                      goto LTM_ERR;
   va_end(args);

   mp_clear(&cur_arg);
LTM_ERR:
   mp_clear(&cur_arg);
   return err;
}

#endif
