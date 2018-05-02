#include <tommath.h>
#ifdef BN_MP_SQRT_C
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
 * Tom St Denis, tomstdenis@gmail.com, http://libtom.org
 */

/* this function is less generic than mp_n_root, simpler and faster */
int mp_sqrt(mp_int *arg, mp_int *ret)
{
   int res, ilog2;
   mp_int t1,t2;

   /* must be positive */
   if (arg->sign == MP_NEG) {
      return MP_VAL;
   }

   /* easy out */
   if (mp_iszero(arg) == MP_YES) {
      mp_zero(ret);
      return MP_OKAY;
   }

   if ((res = mp_init(&t2)) != MP_OKAY) {
      goto E2;
   }

   ilog2 = mp_count_bits(arg) >> 1;
   /* floor(sqrt(n)) == 1 for n in {1,2,3}*/
   if (ilog2 <= 1) {
      mp_set(ret,1);
      goto E1;
   }
   /* First approximation must be bigger than sqrt(arg) to make it converge provably*/
   if ((res = mp_2expt(&t1, ilog2 + 1)) != MP_OKAY) {
      goto E1;
   }
   do {
      if ((res = mp_div(arg,&t1,&t2,NULL)) != MP_OKAY) {
         goto E1;
      }
      if ((res = mp_add(&t1,&t2,&t1)) != MP_OKAY) {
         goto E1;
      }
      if ((res = mp_div_2(&t1,&t1)) != MP_OKAY) {
         goto E1;
      }
   } while (mp_cmp_mag(&t1,&t2) == MP_GT);

   mp_exch(&t1,ret);

E1:
   mp_clear(&t2);
E2:
   mp_clear(&t1);
   return res;
}

#endif

/* $Source$ */
/* $Revision$ */
/* $Date$ */
