#include "tommath_private.h"
#ifdef S_MP_ISPRIME_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



/* deterministic prime test for input < 2^64 */
bool s_mp_isprime_d(mp_digit number)
{
   /* More witnesses (12) than absolutely necessary (7) for testing numbers
      up to 2^64 but all are below 2^16 */
   mp_digit witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

   bool result = false;

   if (number < 2) {
      return false;
   }
   if ((number & 1) == 0) {
      return false;
   }

   /* Instead of small trial divisions  */
   for (int i = 1; i < (int)(sizeof(witnesses)/sizeof(witnesses[0])); i++) {
      if (number == witnesses[i]) {
         return true;
      }
   }

   for (int i = 0; i < (int)(sizeof(witnesses)/sizeof(witnesses[0])) ; i++) {
      result = s_mp_miller_rabin_d(number, witnesses[i]);
      if (!result) {
         return false;
      }
   }
   return true;
}


#endif
