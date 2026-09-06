#include "tommath_private.h"
#ifdef MP_HAMDIST_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

mp_err mp_hamdist(const mp_int *a, const mp_int *b, int *hd)
{
   mp_err err = MP_OKAY;
   mp_int tmp;
   int lena, lenb;
   bool different_length = false;

   lena = (mp_iszero(a)) ? 1 : mp_count_bits(a);
   lenb = (mp_iszero(b)) ? 1 : mp_count_bits(b);

   if (lena != lenb) {
      different_length = true;
   }

   /* Only one of the inputs needs to be tested for zero
     the value of the other one does not matter */
   if (mp_iszero(a)) {
      *hd = mp_popcount(b);
      /* The Hamming distance is normaly between
         two sequences of the same length. To make
         it more flexible it returns MP_OVF (overflow)
         if the input sequences are of different lengths.
         Levenshtein distance is not implemented.
       */
      if (different_length) {
         err = MP_OVF;
      }
      return err;
   }
   if (mp_iszero(b)) {
      *hd = mp_popcount(a);
      if (different_length) {
         err = MP_OVF;
      }
      return err;
   }
   if (!different_length) {
      if (mp_cmp(a,b) == MP_EQ) {
         *hd = 0;
         return err;
      }
   }
   if ((err = mp_init(&tmp)) != MP_OKAY) {
      return err;
   }
   if ((err = mp_xor(a, b, &tmp)) != MP_OKAY)        goto LTM_ERR;
   *hd = mp_popcount(&tmp);
   if (different_length) {
      err = MP_OVF;
   }

LTM_ERR:
   mp_clear(&tmp);
   return err;
}

#endif
