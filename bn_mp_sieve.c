#include "tommath_private.h"
#ifdef BN_MP_SIEVE_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

/* Segmented version of an Eratosthenes sieve */
#define MP_SIEVE_PRIME_NUM_BITS (sizeof(mp_sieve_prime)*CHAR_BIT)

static void s_mp_sieve_setall(mp_single_sieve *bst);
static void s_mp_sieve_clear(mp_single_sieve *bst, mp_sieve_prime n);
static mp_sieve_prime s_mp_sieve_get(mp_single_sieve *bst, mp_sieve_prime n);
static mp_sieve_prime s_mp_sieve_nextset(mp_single_sieve *bst, mp_sieve_prime n);

static mp_sieve_prime s_isqrt(mp_sieve_prime n);

static void s_mp_eratosthenes(mp_single_sieve *bst);
static int  s_mp_eratosthenes_init(mp_sieve_prime n, mp_single_sieve *bst);

static void s_mp_eratosthenes_segment(mp_sieve_prime a, mp_sieve_prime b, mp_single_sieve *base,
                                      mp_single_sieve *segment);
static int s_mp_eratosthenes_segment_init(mp_sieve_prime a, mp_sieve_prime b, mp_single_sieve *base,
      mp_single_sieve *segment);
static void s_mp_eratosthenes_segment_clear(mp_single_sieve *segment, mp_sieve_prime *single_segment_a);

static int s_init_base_sieve(mp_single_sieve *base_sieve);
static int s_init_single_segment_with_start(mp_sieve_prime a, mp_single_sieve *base_sieve,
      mp_single_sieve *single_segment, mp_sieve_prime *single_segment_a);

static mp_sieve_prime s_nextprime(mp_sieve_prime p, mp_single_sieve *bst);


/* Manual memset to avoid the inclusion of string.h */
static void s_mp_sieve_setall(mp_single_sieve *bst)
{
   mp_sieve_prime i, bs_size;
   bs_size = bst->alloc / sizeof(mp_sieve_prime);
   for (i = 0; i < bs_size; i++) {
      (bst)->content[i] = MP_SIEVE_PRIME_MAX;
   }
}

static void s_mp_sieve_clear(mp_single_sieve *bst, mp_sieve_prime n)
{
   ((*((bst)->content+(n/MP_SIEVE_PRIME_NUM_BITS))
     &= ~(1lu<<(n % MP_SIEVE_PRIME_NUM_BITS))));
}

static mp_sieve_prime s_mp_sieve_get(mp_single_sieve *bst, mp_sieve_prime n)
{
   return (((*((bst)->content+(n/MP_SIEVE_PRIME_NUM_BITS))
             & (1lu<<(n % MP_SIEVE_PRIME_NUM_BITS))) != 0));
}

static mp_sieve_prime s_mp_sieve_nextset(mp_single_sieve *bst, mp_sieve_prime n)
{
   while ((n < ((bst)->size)) && (!s_mp_sieve_get(bst, n))) {
      n++;
   }
   return n;
}


/*
 * Integer square root, hardware style
 * Wikipedia calls it "Digit-by-digit calculation"
 * https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Digit-by-digit_calculation
 * This is the base 2 method described at
 * https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Binary_numeral_system_(base_2)
 */
static mp_sieve_prime s_isqrt(mp_sieve_prime n)
{
   mp_sieve_prime s, rem, root;

   if (n < 1uL) {
      return 0uL;
   }
   /* highest power of four <= n */
   s = (mp_sieve_prime)(1uL << (MP_SIEVE_PRIME_NUM_BITS - 2));

   rem = n;
   root = 0uL;
   while (s > n) {
      s >>= 2uL;
   }
   while (s != 0uL) {
      if (rem >= (s | root)) {
         rem -= (s | root);
         root >>= 1uL;
         root |= s;
      } else {
         root >>= 1uL;
      }
      s >>= 2uL;
   }
   return root;
}

/*
 * Set range_a_b to sqrt(MP_SIEVE_PRIME_MAX)
 * TODO: Make it const or put it in bncore.c because it is said to be faster
 * if the size of range_a_b fits into the L2-cache.
 * Not much difference on the author's machine for 32 bit but quite
 * a large one for 64 bit and large limits. YMMV, as always.
 * Please be aware that range_a_b is in bits, not bytes and memory
 * allocation rounds up and adds CHAR_BIT*sizeof(mp_sieve_prime) bits.
 */
#ifndef MP_SIEVE_RANGE_A_B
#if ( (defined MP_64BIT) && (defined MP_SIEVE_USE_LARGE_SIEVE) )
#define MP_SIEVE_RANGE_A_B  0x400000uL
#else
#define MP_SIEVE_RANGE_A_B ((mp_sieve_prime) MP_SIEVE_PRIME_MAX_SQRT)
#endif
#endif
#define MP_SIEVE_BASE_SIEVE_SIZE  ((mp_sieve_prime)MP_SIEVE_PRIME_MAX_SQRT)

/* TODO: Some redundant code below, needs a clean-up */


/*
 * Initiate a sieve that stores the odd numbers only:
 * allocate memory, set actual size and allocated size and fill it
 */
static int s_mp_eratosthenes_init(mp_sieve_prime n, mp_single_sieve *bst)
{
   n = (n - 1uL) / 2uL;
   bst->content =
      MP_MALLOC(((size_t)n + sizeof(mp_sieve_prime)) / sizeof(mp_sieve_prime) + sizeof(mp_sieve_prime));
   if (bst->content == NULL) {
      return MP_MEM;
   }

   bst->alloc =
      (n + sizeof(mp_sieve_prime)) / sizeof(mp_sieve_prime) + sizeof(mp_sieve_prime);
   bst->size = (2uL * n) + 1uL;

   s_mp_eratosthenes(bst);
   return MP_OKAY;
}

/*
 * Simple Eratosthenes' sieve, starting at zero
 * Keeping odd primes only as the single optimization
 */
static void s_mp_eratosthenes(mp_single_sieve *bst)
{
   mp_sieve_prime n, k, r, j;

   n = bst->size;
   r = s_isqrt(n);
   s_mp_sieve_setall(bst);
   for (k = 1uL; k < ((r - 1uL) / 2uL); k += 1uL) {
      if (s_mp_sieve_get(bst, k)) {
         for (j = k * 2uL * (k + 1uL); j < (n - 1uL) / 2uL; j += 2uL * k + 1uL) {
            s_mp_sieve_clear(bst, j);
         }
      }
   }
}

/*
 * For a sieve that has only the odd numbers stored.
 * It returns next prime even if p is prime in contrast to
 * e.g.: Pari/GP's nextprime() function
 * Used internally for the segment to get the primes from
 * the base sieve.
 */
static mp_sieve_prime s_nextprime(mp_sieve_prime p, mp_single_sieve *bst)
{
   mp_sieve_prime ret;
   if (p <= 1uL) {
      return 2uL;
   }
   ret = s_mp_sieve_nextset(bst, ((p - 1uL) / 2uL) + 1uL);
   return (2uL * ret) + 1uL;
}

/*
 * Init sieve "segment" of the range [a,b]:
 * allocate memory, set actual size and allocated size and fill it from sieve "base"
 * TODO: merge with s_mp_eratosthenes_init()
 */
static int s_mp_eratosthenes_segment_init(mp_sieve_prime a, mp_sieve_prime b,
      mp_single_sieve *segment, mp_single_sieve *base)
{
   mp_sieve_prime n;

   n = b - a;

   segment->content =
      MP_MALLOC(((size_t)n + sizeof(mp_sieve_prime)) / sizeof(mp_sieve_prime) +
                sizeof(mp_sieve_prime));
   if (segment->content == NULL) {
      return MP_MEM;
   }
   segment->alloc =
      (n + sizeof(mp_sieve_prime)) / sizeof(mp_sieve_prime) +
      sizeof(mp_sieve_prime);
   segment->size = n;

   s_mp_eratosthenes_segment(a, b, base, segment);
   return MP_OKAY;
}

/*
 * Clear sieve "segment"
 * free memory and reset "single_segment_a"
 */
static void s_mp_eratosthenes_segment_clear(mp_single_sieve *segment,
      mp_sieve_prime *single_segment_a)
{
   if (segment->content != NULL) {
      MP_FREE(segment->content, ((size_t)segment->size + sizeof(mp_sieve_prime)) /
              sizeof(mp_sieve_prime) + sizeof(mp_sieve_prime));
   }
   *single_segment_a = 0uL;

}

/*
 * TODO: Not memory optimized, it stores the even numbers, too.
 * Should be ok if the segments are small but needs to
 * be done at some time in the near future.
 *
 * Fill sieve "segment" of the range [a,b] from the basic sieve "base"
 */
static void s_mp_eratosthenes_segment(mp_sieve_prime a, mp_sieve_prime b,
                                      mp_single_sieve *base, mp_single_sieve *segment)
{
   mp_sieve_prime r, j, i;
   mp_sieve_prime p;
   r = s_isqrt(b);
   s_mp_sieve_setall(segment);
   p = 0uL;
   for (i = 0uL; p <= r; i++) {
      p = s_nextprime(p, base);
      j = p * p;
      if (j < a) {
         j = ((a + p - 1uL) / p) * p;
      }
      for (; j <= b; j += p) {
         /* j+=p can overflow */
         if (j >= a) {
            s_mp_sieve_clear(segment, j - a);
         } else {
            break;
         }
      }
   }
}

/* Build a basic sieve with the largest reasonable size */
static int s_init_base_sieve(mp_single_sieve *base_sieve)
{
   mp_sieve_prime n;
   /*
    * That is a risky idea. If range_a_b is smaller than sqrt(n_max)
    * the whole thing stops working!
    */
   /* n = range_a_b; */
   n = MP_SIEVE_BASE_SIEVE_SIZE;
   return s_mp_eratosthenes_init(n, base_sieve);
}

/*
 * Build a segment sieve with the largest reasonable size. "a" is the start of
 * the sieve Size is MIN(range_a_b,MP_SIEVE_PRIME_MAX-a)
 */
static int s_init_single_segment_with_start(mp_sieve_prime a,
      mp_single_sieve *base_sieve,
      mp_single_sieve *single_segment,
      mp_sieve_prime *single_segment_a)
{
   mp_sieve_prime b;
   int e = MP_OKAY;

   /* last segment might not fit, depending on size of range_a_b */
   if (a > (MP_SIEVE_PRIME_MAX - MP_SIEVE_RANGE_A_B)) {
      b = (mp_sieve_prime)MP_SIEVE_PRIME_MAX - 1uL;
   } else {
      b = a + (mp_sieve_prime)MP_SIEVE_RANGE_A_B;
   }
   if ((e =
           s_mp_eratosthenes_segment_init(a, b, single_segment,
                                          base_sieve)) != MP_OKAY) {
      return e;
   }
   *single_segment_a = a;
   return e;
}

/*
 * Sets "result" to one if n is prime or zero respectively.
 * Also sets "result" to zero in case of error.
 * Worst case runtime is: building a base sieve and a segment and
 * search the segment
 */
mp_err mp_is_small_prime(mp_sieve_prime n, mp_sieve_prime *result, mp_sieve *sieve)
{
   int e = MP_OKAY;
   mp_sieve_prime a = 0uL, b = 0uL;

   if (n < 2uL) {
      *result = 0uL;
      return e;
   }

   if ((n & 1uL) == 0uL) {
      *result = 0uL;
      return e;
   }

   /* neither of 2^16-x, 2^32-x, or 2^64-x are prime for 0<=x<=4 */
   if (n >= (mp_sieve_prime)(MP_SIEVE_PRIME_MAX - 3)) {
      *result = 0uL;
      return e;
   }

   if (sieve->base.content == NULL) {
      if ((e = s_init_base_sieve(&(sieve->base))) != MP_OKAY) {
         *result = 0uL;
         return e;
      }
   }

   /* No need to generate a segment if n is in the base sieve */
   if (n < MP_SIEVE_BASE_SIEVE_SIZE) {
      /* might have been a small sieve, so check size of sieve first */
      if (n < sieve->base.size) {
         *result = s_mp_sieve_get(&(sieve->base), (n - 1uL) / 2uL);
         return e;
      }
   }

   /* no further shortcuts to apply, build and search a segment */

   /*
    * TODO: if base_sieve is < sqrt(MP_SIEVE_PRIME_MAX) it is not possible to get
    * all primes <= MP_SIEVE_PRIME_MAX so add check if n > size(base_sieve)^2 and either
    * a. make size(base_sieve) = sqrt(MP_SIEVE_PRIME_MAX)
    * b. make size(base_sieve) = 2 * size(base_sieve)
    * with 2 * size(base_sieve) <= sqrt(MP_SIEVE_PRIME_MAX)
    * c. give up and return MP_VAL
    */
   /* we have a segment and may be able to use it */
   if (sieve->segment.content != NULL) {
      a = sieve->single_segment_a;
      /* last segment may not fit into range_a_b */
      if (a > (MP_SIEVE_PRIME_MAX - MP_SIEVE_RANGE_A_B)) {
         b = (mp_sieve_prime)(MP_SIEVE_PRIME_MAX - 1);
      } else {
         b = a + MP_SIEVE_RANGE_A_B;
      }
      /* check if n is inside the bounds of the segment */
      if (n >= a && n <= b) {
         *result = s_mp_sieve_get(&(sieve->segment), n - a);
         return e;
      }
      /* get a clean slate */
      else {
         s_mp_eratosthenes_segment_clear(&(sieve->segment), &(sieve->single_segment_a));
      }
   }

   /*
    * A bit of heuristics ( "heuristics" is a more pretentious word for the
    * commonly known expression "wild guess")
    * Based on the vague idea of the assumption that most sieves get used for
    * sequential series of primes or a single test here and there, but not
    * for massive amounts of random requests.
    */
   if (n > a) {
      if (n > (MP_SIEVE_PRIME_MAX - MP_SIEVE_RANGE_A_B)) {
         a = (mp_sieve_prime)(MP_SIEVE_PRIME_MAX - MP_SIEVE_RANGE_A_B);
      } else {
         a = n;
      }
   } else {
      if (n < MP_SIEVE_RANGE_A_B) {
         a = (mp_sieve_prime)(MP_SIEVE_BASE_SIEVE_SIZE - 1);
      } else {
         /* TODO: there should be no need for an overlap, check */
         a = n - (mp_sieve_prime)(MP_SIEVE_RANGE_A_B - 2);
      }
   }
   if ((e = s_init_single_segment_with_start(a,
            &(sieve->base), &(sieve->segment), &(sieve->single_segment_a))) != MP_OKAY) {
      *result = 0uL;
      return e;
   }
   /* finally, check for primality */
   *result = s_mp_sieve_get(&(sieve->segment), n - a);
   return e;
}
#endif
