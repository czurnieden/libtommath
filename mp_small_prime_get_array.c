#include "tommath_private.h"
#ifdef MP_SMALL_PRIME_GET_ARRAY_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */




mp_err mp_small_prime_get_array(const mp_digit start, const mp_digit end, int *count, ERAT_UINT *last_prime,
                                ERAT_UINT **prime_array, mp_erat_sieve *sieve)
{
   mp_err err = MP_OKAY;
   int counter = 1;
   size_t allocated = 0u, i = 0u;
   bool want_count = false, wantlastprime = false, wantarray_filled = false;

   ERAT_UINT p, p_last;
   ERAT_UINT *tmp = NULL;

   if (end < start) {
      return MP_VAL;
   }
   /* The size of ERAT_UINT is fixed, the size of mp_digit is not */
#if (MP_DIGIT_BIT > 16)
   if ((start > (mp_digit)ERAT_UINT_MAX) || (end > (mp_digit)ERAT_UINT_MAX)) {
      return MP_VAL;
   }
#endif
   if (end < 2) {
      return MP_VAL;
   }

   if (count != NULL) {
      want_count = true;
   }
   if (last_prime != NULL) {
      wantlastprime = true;
   }
   if (prime_array != NULL) {
      wantarray_filled = true;
      /* Or free() it? Or add to it? */
      if (*prime_array != NULL) {
         return MP_VAL;
      }
   }

   /* User cannot want nothing, that's an absolute no-go!
      Where would it head to if we endorse such odd behaviour? */
   if (!want_count && !wantlastprime && !prime_array) {
      return MP_VAL;
   }


   /* A bit expensive, but if the caller insists on it ... */
   if (!want_count && wantlastprime && !prime_array) {
      /*
          That will be wrong if the difference "end"-"start"
          is smaller than the gap to the prime preceding "end".
       */
      if ((err = mp_small_prime_sieve_prec_prime((ERAT_UINT)end - 1u, &p, sieve)) != MP_OKAY) {
         return err;
      }
      if ((mp_digit)p < start) {
         return MP_VAL;
      }
      *last_prime = p;
      return err;
   }



   /* "start" might be a prime  */
   if ((err = mp_small_prime_sieve_next_prime(start, &p, sieve)) != MP_OKAY) {
      return err;
   }
   if (p > (ERAT_UINT)end) {
      return MP_VAL;
   }
   if (wantarray_filled) {
      *prime_array = MP_MALLOC((size_t)MP_SMALL_PRIME_COUNT_REALLOC_STEP * sizeof(ERAT_UINT));
      if (*prime_array == NULL) {
         return MP_MEM;
      }
      allocated = (size_t)MP_SMALL_PRIME_COUNT_REALLOC_STEP;

      (*prime_array)[i] = p;
      i++;
   }
   while (true) {
      p_last = p;
      if ((err = mp_small_prime_sieve_next_prime(p + 1u, &p, sieve)) != MP_OKAY) {
         return err;
      }
      if (p > end) {
         if (wantlastprime) {
            *last_prime = (uint64_t)p_last;
         }
         break;
      }
      if (wantarray_filled) {
         (*prime_array)[i] = p;
         i++;
      }
      /* TODO: Counts primes between start and end, don't forget to take a note in the docs */
      counter++;
      if (wantarray_filled) {
         if (i >= allocated) {
            tmp = MP_REALLOC(*prime_array, allocated * sizeof(ERAT_UINT),
                             (allocated + (size_t)MP_SMALL_PRIME_COUNT_REALLOC_STEP) * sizeof(ERAT_UINT));
            if (tmp == NULL) {
               MP_FREE(*prime_array, allocated * sizeof(ERAT_UINT));
               return MP_MEM;
            }
            allocated += (size_t)MP_SMALL_PRIME_COUNT_REALLOC_STEP;
            *prime_array = tmp;
         }
      }
   }
   if (wantarray_filled) {
      tmp = MP_REALLOC(*prime_array, (size_t)counter * sizeof(ERAT_UINT), (size_t)counter * sizeof(ERAT_UINT));
      if (tmp == NULL) {
         MP_FREE(*prime_array, (size_t)counter * sizeof(ERAT_UINT));
         return MP_MEM;
      }
      *prime_array = tmp;
   }
   if (want_count) {
      *count = counter;
   }
   return err;
}



#endif
