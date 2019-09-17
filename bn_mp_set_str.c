#include "tommath_private.h"
#ifdef BN_MP_SET_STR_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

/*
   Simple D&C algorithm to speed up number conversion. ( O(log(n) * M(n))) )
   Stand-alone, can replace mp_read_radix.
*/

/* TODO: change to "size_t" when the type of "mpint.size" changes to "size_t" */
static int s_mp_strlen(const char *s)
{
   const char *p;
   p = s;
   while (*p != '\0') {
      p++;
   }
   return (int)(p - s);
}

/* ASCII only */
#if  ('\n' == 0x0A) && (' ' == 0x20) && ('0' == 0x30) && ('A' == 0x41) && ('a' == 0x61) \
     && ('!' == 0x21)
#define MP_IS_ASCII
#else
#error "ASCII only for now, sorry."
#endif

/* no ASCII characters "lower" than 0x28 which is a "(", an open parentheses */
#ifndef MP_ISUSELESS
#define MP_ISUSELESS(c)  ( (c) < 0x28 )
#endif

#ifndef MP_TOUPPER
#define MP_TOUPPER(c) ((((c) >= 'a') && ((c) <= 'z')) ? (((c) + 'A') - 'a') : (c))
#endif

/*
   TODO:
   This can be replaced with
      mp_err mp_read_radix(mp_int *a, const char *str, int radix) {
         return s_mp_read_radix_chunk(a, str, 0, s_mp_strlen(str), radix);
      }
   if The Powers That Be can be persuaded.
*/

static mp_err s_mp_set_str_chunk(mp_int *a, const char *str, const int start, const int end, int base)
{
   mp_err err;
   char *_s;
   char ch;
   unsigned int pos;
   int i, y;

   /* checks are done by caller */

   _s = (char *)(str + start);
   for (i = start; i < end; i++) {
      if (MP_ISUSELESS((int)(*_s))) {
         err = MP_VAL;
         goto LBL_ERR;
      }
      /* description in bn_mp_read_radix.c */
      ch = (base <= 36) ? (char)MP_TOUPPER((int)*_s) : *_s;
      pos = (unsigned)(ch - '(');
      y = (int)mp_s_rmap_reverse[pos];
      if ((y == 0xff) || (y >= base)) {
         err = MP_VAL;
         goto LBL_ERR;
      }
      if ((err = mp_mul_d(a, (mp_digit)base, a)) != MP_OKAY)      goto LBL_ERR;
      if ((err = mp_add_d(a, (mp_digit)y, a)) != MP_OKAY)         goto LBL_ERR;
      _s++;
   }

   return MP_OKAY;
LBL_ERR:
   mp_zero(a);
   return err;
}


/* Needs to get computed based on radix (table?), this value is for base 10 */
static int SET_STR_CUTOFF = 100;

static mp_err s_mp_set_str_intern(mp_int *a, const char *string, const int start, const int end, int base)
{
   int len, mid;
   mp_int A, B, m;
   mp_err e = MP_OKAY;

   len = end - start;

   if (len < SET_STR_CUTOFF) {
      return s_mp_set_str_chunk(a, string, start, end, base);
   }

   mid = len / 2;

   if ((e = mp_init_set(&m, base)) != MP_OKAY) {
      return e;
   }
   if ((e = mp_init_multi(&A, &B, NULL)) != MP_OKAY) {
      mp_clear(&m);
      return e;
   }
   /*
       The chunk of "string we" are working at is between "start" and "end" where
       both "start" and "end" are values relative to the orginal string so with the
       very first call of this function "start" is the start of the original "string"
       (or 0 (zero)) and "end" is the end of the original "string" without the final
       '\0' (strlen(string)).

       The divides are quite simple, just cut the length in half, so we get
       (with {A, B} the partial binary numbers) the new positions

          A) new_start = start,           new_end = start + mid
          B) new_start = start + mid + 1, new_end = end

       with mid = floor((end - start)/2)
   */

   if ((e = s_mp_set_str_intern(&A, string, start, start + mid + 1, base)) != MP_OKAY) {
      goto LBL_ERR;
   }
   if ((e = s_mp_set_str_intern(&B, string, start + mid + 1, end, base)) != MP_OKAY) {
      goto LBL_ERR;
   }
   /* TODO: caching? */
   /* TODO: replace with shifts for bases 2, 4, 8, 16, 64 or change s_mp_set_str_chunk
             for the cost of a larger binary. */
   if ((e = mp_expt_u32(&m, (uint32_t)((len - mid) - 1), &m)) != MP_OKAY)     goto LBL_ERR;
   if ((e = mp_mul(&A, &m, &A)) != MP_OKAY)                             goto LBL_ERR;
   if ((e = mp_add(&A, &B, a)) != MP_OKAY)                              goto LBL_ERR;

LBL_ERR:
   mp_clear_multi(&A, &B, &m, NULL);
   return e;
}

mp_err mp_set_str(mp_int *a, const char *string, int base)
{
   mp_sign sign;
   mp_err err;
   int len;

   if ((base < 2) || (base > 64)) {
      return MP_VAL;
   }

   mp_zero(a);

   len = s_mp_strlen(string);

   /* TODO: trim string? */
   if (*string == '-') {
      string++;
      len--;
      sign = MP_NEG;
   } else {
      sign = MP_ZPOS;
   }
   /* A small shortcut */
   if (len <= SET_STR_CUTOFF) {
      if ((err = s_mp_set_str_chunk(a, string, 0, len, base)) != MP_OKAY) {
         return err;
      }
   } else {
      if ((err = s_mp_set_str_intern(a, string, 0, len, base)) != MP_OKAY) {
         return err;
      }
   }
   if (!MP_IS_ZERO(a)) {
      a->sign = sign;
   }

   return MP_OKAY;
}

#endif
