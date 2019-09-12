#include "tommath_private.h"
#ifdef BN_MP_SET_STR_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



/*
   Simple D&C algorithm to speed up number conversion
   For comments see mp_get_str
*/

#include <string.h>

static int SET_STR_CUTOFF = 20;
static mp_err s_mp_set_str_intern(mp_int *a, const char *string, int base)
{
   int len, len_low, len_high;
   char *s_low, *s_high;
   mp_int A, B, m;
   mp_err e = MP_OKAY;

   len = strlen(string);
   len_low = len / 2;
   len_high = len - len_low;

   if (len_low < SET_STR_CUTOFF) {
      return mp_read_radix(a, string, base);
   }
   if ((e = mp_init_set(&m, base)) != MP_OKAY) {
      return e;
   }
   if ((e = mp_init_multi(&A, &B, NULL)) != MP_OKAY) {
      mp_clear(&m);
      return e;
   }
   /*
      Yes, but do _you_ want to do all of the necessary pointer
      juggling for cheesy little prototype?
   */
   s_high = malloc(sizeof(char) * (len_high + 1));
   if (s_high == NULL) {
      e = MP_MEM;
      goto LBL_ERR;
   }
   memcpy(s_high, string + len_low, len_high);
   s_high[len_high] = '\0';
   if ((e = s_mp_set_str_intern(&A, s_high, base)) != MP_OKAY) {
      free(s_high);
      goto LBL_ERR;
   }
   free(s_high);

   s_low = malloc(sizeof(char) * (len_low + 1));
   if (s_low == NULL) {
      e = MP_MEM;
      goto LBL_ERR;
   }
   memcpy(s_low, string, len_low);
   s_low[len_low] = '\0';
   if ((e = s_mp_set_str_intern(&B, s_low, base)) != MP_OKAY) {
      free(s_low);
      goto LBL_ERR;
   }
   free(s_low);
   /* TODO: caching */
   if ((e = mp_expt_u32(&m, (uint32_t)len_high, &m)) != MP_OKAY)        goto LBL_ERR;
   if ((e = mp_mul(&B, &m, &B)) != MP_OKAY)                             goto LBL_ERR;
   if ((e = mp_add(&A, &B, a)) != MP_OKAY)                              goto LBL_ERR;

LBL_ERR:
   mp_clear_multi(&A, &B, &m, NULL);
   return e;
}

mp_err mp_set_str(mp_int *a, const char *string, int base)
{
   mp_sign sign;
   mp_err err;
   /* TODO: trim string? */
   if (*string == '-') {
      string++;
      sign = MP_NEG;
   } else {
      sign = MP_ZPOS;
   }
   if (*string == '+') {
      string++;
   }
   if ((err = s_mp_set_str_intern(a, string, base)) != MP_OKAY) {
      a->sign = sign;
      return err;
   }
   a->sign = sign;
   return MP_OKAY;
}



#endif
