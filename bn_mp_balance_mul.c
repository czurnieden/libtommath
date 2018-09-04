#include <tommath_private.h>
#ifdef BN_MP_BALANCE_MUL_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis
 *
 * LibTomMath is a library that provides multiple-precision
 * integer arithmetic as well as number theoretic functionality.
 *
 * The library was designed directly after the MPI library by
 * Michael Fromberger but has been written from scratch with
 * additional optimizations in place.
 *
 * The library is free for all purposes without any express
 * guarantee it works.
 *
 * Tom St Denis, tstdenis82@gmail.com, http://libtom.org
 */

/*
 *  Balancing multiplication for Toom-Cook algorithms and alike.
 *  Implemented as a simple multi-digit by single-digit multiplication
 *  with the smaller big integer as the single-digit
 */
int mp_balance_mul(const mp_int *a, const mp_int *b, mp_int *c)
{
   int e, count, nblocks, i, j, bsize;
   mp_int a0, tmp, r;
   const mp_int *A, *B;

   nblocks = MAX(a->used, b->used) / MIN(a->used, b->used);
   bsize = MIN(a->used, b->used) ;
   e = MP_OKAY;

   if ((e = mp_init_size(&a0, bsize + 1)) != MP_OKAY) {
      return e;
   }
   if ((e = mp_init_multi(&tmp, &r, NULL)) != MP_OKAY) {
      mp_clear(&a0);
      return e;
   }

   /* Make sure that A is the larger one */
   if (a->used < b->used) {
      B = a;
      A = b;
   } else {
      A = a;
      B = b;
   }

   /*
      These blocks can get very large and hence use a lot of
      heap-memory.
   */
   for (i = 0, j=0; i < nblocks; i++) {
      a0.used = 0;
      for (count = 0; count < bsize; count++) {
         a0.dp[count] = A->dp[ j++ ];
         a0.used++;
      }
      if ((e = mp_mul(&a0, B, &tmp)) != MP_OKAY) {
         goto ERR;
      }
      if ((e = mp_lshd(&tmp, bsize * i)) != MP_OKAY) {
         goto ERR;
      }
      /* Just add to output. No carry needed */
      if ((e = mp_add(&r, &tmp, &r)) != MP_OKAY) {
         goto ERR;
      }
   }
   /* The left-overs; there are almost always left-overs */
   if (j < A->used) {
      a0.used = 0;
      for (count = 0; j < A->used; count++) {
         a0.dp[count] = A->dp[ j++ ];
         a0.used++;
      }
      if ((e = mp_mul(&a0, B, &tmp)) != MP_OKAY) {
         goto ERR;
      }
      if ((e = mp_lshd(&tmp, bsize * i)) != MP_OKAY) {
         goto ERR;
      }
      if ((e = mp_add(&r, &tmp, &r)) != MP_OKAY) {
         goto ERR;
      }
   }

   mp_exch(&r,c);
ERR:
   mp_clear_multi(&a0, &tmp, &r,NULL);
   return e;
}
#endif
