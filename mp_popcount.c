#include "tommath_private.h"
#ifdef MP_POPCOUNT_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */


/* MSVC offers a built-in of a similar-looking name but it does not work the same.
   But if you have the knowledge and the patience...

   The macro __popcount AMD is in intrin.h with following prototypes:

      unsigned short __popcnt16(unsigned short value);
      unsigned int __popcnt(unsigned int value);
      unsigned __int64 __popcnt64(unsigned __int64 value);

   The macro _mm_popcnt_?? (Intel specific) is in intrin.h, too, with the following prototypes:

      int _mm_popcnt_u32 (unsigned int a)
      __int64 _mm_popcnt_u64 (unsigned __int64 a)

   Those two are binary compatible with the __popcnt?? functions.

   It seems as if MSVC uses the assembler operation without checking if
   the architecture supports it. Use cpuid to check it yourself:

      void __cpuid(int cpuInfo[4], int function_id);

   Snippet (for both versions, AMD and Intel, no need to check for SSE4-2 support):

      #include <intrin.h>
      int check_popcnt(void) {
         int cpu_info[4], ecx;
         __cpuid(cpu_info, 0x1);
         ecx = cpu_info[2];
         return (ecx & (1 << 23)) != 0;
      }

 */
#if defined(__GNUC__) || defined(__clang__)
#ifdef MP_16BIT
#define MP_CAST_TYPE  (unsigned int)
#define bltn_popcount(x)  __builtin_popcount((x))
#elif ((defined MP_32BIT) || (defined MP_31BIT) || (defined MP_28BIT) )
#define MP_CAST_TYPE  (unsigned long)
#define bltn_popcount(x)  __builtin_popcountl((x))
#elif (defined MP_64BIT)
#define MP_CAST_TYPE  (unsigned long long)
#define bltn_popcount(x)  __builtin_popcountll((x))
#elif
#error "No digit-size (MP_??BIT) defined, please get in contact"
#endif

/* Otherwise use the common divide&conquer method */
#else

#ifdef MP_16BIT
#define MP_ALL_FIVERS 0x5555u
#define MP_ALL_THREES 0x3333u
#define MP_ALL_OUGHFS 0x0F0Fu
#define MP_ALL_OHONES 0x0101u
#elif ((defined MP_32BIT) || (defined MP_31BIT) || (defined MP_28BIT) )
#define MP_ALL_FIVERS 0x55555555ul
#define MP_ALL_THREES 0x33333333ul
#define MP_ALL_OUGHFS 0x0F0F0F0Ful
#define MP_ALL_OHONES 0x01010101ul
#elif (defined MP_64BIT)
#define MP_ALL_FIVERS 0x5555555555555555ull
#define MP_ALL_THREES 0x3333333333333333ull
#define MP_ALL_OUGHFS 0x0F0F0F0F0F0F0F0Full
#define MP_ALL_OHONES 0x0101010101010101ull
#elif
#error "No digit-size (MP_??BIT) defined, please get in contact"
#endif

#endif
/*
   We can do this trick with the whole bigint at once, but only if we need it often,
   the bigints have always the same limb-sizes and only if we cache the magik numbers.
   It is not only the generation of the constants, we need to haul the whole bigint
   to and fro for all of the computations. Not much individually, but it adds.
   It is still O(log n) but the constants gets large.

   Doing it with native small integers should be faster on many archtitectures.

   BUT: YMMV as always.
*/


int mp_popcount(const mp_int *a)
{
   int total_set_bits = 0, i;
#if !(defined(__GNUC__) || defined(__clang__))
   mp_digit x0, x1, x2, x3, allx;
#endif
   if (mp_iszero(a)) {
      return 0;
   }

   /* Does loop-unrolling make sense in general? With modern compilers, too? */
   for (i = 0; i <= a->used - 4; i += 4) {
#if defined(__GNUC__) || defined(__clang__)
      total_set_bits += bltn_popcount(MP_CAST_TYPE a->dp[i]);
      total_set_bits += bltn_popcount(MP_CAST_TYPE a->dp[i + 1]);
      total_set_bits += bltn_popcount(MP_CAST_TYPE a->dp[i + 2]);
      total_set_bits += bltn_popcount(MP_CAST_TYPE a->dp[i + 3]);
#else
      x0 = a->dp[i];
      x1 = a->dp[i+1];
      x2 = a->dp[i+2];
      x3 = a->dp[i+3];

      x0 -= (x0 >> 1) & MP_ALL_FIVERS;
      x1 -= (x1 >> 1) & MP_ALL_FIVERS;
      x2 -= (x2 >> 1) & MP_ALL_FIVERS;
      x3 -= (x3 >> 1) & MP_ALL_FIVERS;

      x0 = (x0 & MP_ALL_THREES) + ((x0 >> 2) & MP_ALL_THREES);
      x1 = (x1 & MP_ALL_THREES) + ((x1 >> 2) & MP_ALL_THREES);
      x2 = (x2 & MP_ALL_THREES) + ((x2 >> 2) & MP_ALL_THREES);
      x3 = (x3 & MP_ALL_THREES) + ((x3 >> 2) & MP_ALL_THREES);

      x0 = (x0 + (x0 >> 4)) & MP_ALL_OUGHFS;
      x1 = (x1 + (x1 >> 4)) & MP_ALL_OUGHFS;
      x2 = (x2 + (x2 >> 4)) & MP_ALL_OUGHFS;
      x3 = (x3 + (x3 >> 4)) & MP_ALL_OUGHFS;

      allx = x0 + x1 + x2 + x3;

      total_set_bits += (int)((allx * MP_ALL_OHONES) >> 56);
#endif
   }
   for (; i < a->used; i++) {
#if defined(__GNUC__) || defined(__clang__)
      total_set_bits += bltn_popcount((unsigned long long)a->dp[i]);
#else
      uint64_t x = a->dp[i];
      x -= (x >> 1) & MP_ALL_FIVERS;
      x = (x & MP_ALL_THREES) + ((x >> 2) & MP_ALL_THREES);
      x = (x + (x >> 4)) & MP_ALL_OUGHFS;
      total_set_bits += (int)((x * MP_ALL_OHONES) >> 56);
#endif
   }
   return total_set_bits;
}

#endif
