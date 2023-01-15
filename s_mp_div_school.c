#include "tommath_private.h"
#ifdef S_MP_DIV_SCHOOL_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

/* integer signed division.
 * c*b + d == a [e.g. a/b, c=quotient, d=remainder]
 * HAC pp.598 Algorithm 14.20
 *
 * Note that the description in HAC is horribly
 * incomplete.  For example, it doesn't consider
 * the case where digits are removed from 'x' in
 * the inner loop.  It also doesn't consider the
 * case that y has fewer than three digits, etc..
 *
 * The overall algorithm is as described as
 * 14.20 from HAC but fixed to treat these cases.
*/

#ifdef USE_KNUTH
mp_err s_mp_div_school(const mp_int *u, const mp_int *v, mp_int *Q, mp_int *R)
{
   const mp_word base = MP_DIGIT_MAX + 1;
   mp_int un, vn, q;
   mp_word qhat, rhat, tmpw;
   mp_digit borrow, tmpd;
   bool neg;

   int s;
   int i, j;
   int m, n;

   mp_err err = MP_OKAY;

   neg = (u->sign != v->sign);

   m = u->used;
   n = v->used;

   if ((err = mp_init_size(&q, m + 2)) != MP_OKAY) {
      return err;
   }
   q.used = m + 2;

   /* We need a bit of leeway here. */
   if ((err = mp_init_size(&un, m + 2)) != MP_OKAY)        goto LTM_ERR_0;
   if ((err = mp_init_size(&vn, n + 2)) != MP_OKAY)        goto LTM_ERR;

   /*
      D1. [Normalize]
         Set D to (B - 1) / V[n - 1]
    */
   s = mp_count_bits(v) % MP_DIGIT_BIT;
   if (s != 0) {
      s = MP_DIGIT_BIT - s;
   }
   if ((err = mp_mul_2d(u, s, &un)) != MP_OKAY)             goto LTM_ERR;
   if ((err = mp_mul_2d(v, s, &vn)) != MP_OKAY)             goto LTM_ERR;
   un.dp[un.used] = 0u;
   vn.dp[vn.used] = 0u;

   /*
      D2. [Initialize j]
          Set the loop counter j to m.
    */
   for (j = m - n; j >= 0; j--) {
      /*
        D3. [Calculate Qhat]
            Set Qhat to (U[n+j] x B + U[n-1+j]) / V[n-1];
            Set Rhat to (U[n+j] x B + U[n-1+j]) % V[n-1];
            Test if Qhat equals B or Qhat * V[n-2] is greater than Rhat * B + U[n-2+j];
            If yes, then decrease Qhat by 1, increase Rhat by V[n-1], and repeat this test while R is less than B.

            Short: the first approximation of the quotient/remainder is made by dividing the first two digits
            of the current numerator by the first digit of the current denominator. Use a third digit from the
            current numerator and second from denominator to refine that approximation
      */
      tmpw = ((mp_word)un.dp[j + n] << MP_DIGIT_BIT) | ((mp_word)un.dp[j + n - 1]);
      qhat = tmpw / (mp_word)vn.dp[n - 1];
      rhat = tmpw % (mp_word)vn.dp[n - 1];

      /* Otherwise "n - 2 < 0"  */
      if (vn.used > 1) {
         for (;;) {
            if ((rhat < base) && (qhat * vn.dp[n - 2]) > (base * rhat + un.dp[j + n - 2])) {
               qhat--;
               rhat = rhat + vn.dp[n - 1];
            } else {
               break;
            }
         }
      }

      /*
        D4. [Multiply and subtract]
           Replace (U[n+j]U[n-1+j]...U[j]) by (U[n+j]U[n-1+j]...U[j]) - Qhat * (V[n-1]...V[1]V[0]).
           (The "digits" (U[n+j}...U[j]) should be kept positive; if the result of this step is actually negative,
           (U[n+j]...U[j]) should be left as the true value plus Bn+1, namely as the B's complement of the true value,
           and a borrow to the left should be remembered.)

           That is easy if the digits (limbs) are of a signed type. We have an unsigned type here and have to go to some
           length to circumnavigate that. Also: the two arrays we are iterating over are of different length but the
           difference is only one, so either add another step after the loop or check inside the loop.
      */

      borrow = 0;
      for (i = 0; i < n; i++) {
         /* Same as above without the temporary variable "p". Not needed here. */
         tmpw = (mp_word)un.dp[i + j] - qhat * (mp_word)vn.dp[i] - (mp_word)borrow;
         /* Produce a complement from the high part. This time from an unsigned integer where
            an unary minus is sufficient. */
         borrow = -((mp_digit)(tmpw >> MP_DIGIT_BIT));
         un.dp[i + j] = (mp_digit)(tmpw & MP_MASK);
      }
      /* Last round with the non-existing vn.dp[n] set to zero */
      tmpw = (mp_word)un.dp[i + j] - (mp_word)borrow;
      borrow = -((mp_digit)(tmpw >> MP_DIGIT_BIT));
      un.dp[i + j] = (mp_digit)(tmpw & MP_MASK);

      /*
         D5. [Test remainder]
            Set Q[j] to Qhat;
            If the result of step D4 was negative, i.e. the subtraction needed a borrow, then proceed with step D6;
            otherwise proceed with step D7.
      */
      q.dp[j] = (mp_digit)(qhat & MP_MASK);

      /*
         D6. [Add back]
             Decrease Q[j] by 1 and add (0V[n-1]...V[1]V[0]) to (U[n+j]U[n-1+j]...U[1+j]U[j]).
             (A carry will occur to the left of U[n+j], and it should be ignored since it cancels with the borrow that occurred in step D4.)
       */

      if (borrow != 0) {
         /* mp_digit is always smaller than the underlying type, so UTYPE_MAX - 1 will get us all bits set
            which is a bit too much, obviously. Error came up after 2.3 mio iterations of random tests. */
         if (q.dp[j] == 0x0) {
            q.dp[j] = MP_MASK;
         } else {
            q.dp[j] = q.dp[j] - 1;
         }

         /* Variable reuse rarely raises readability. */
         borrow = 0;
         for (i = 0; i <= n; i++) {
            tmpw = (mp_word)vn.dp[i] + (mp_word)un.dp[i + j] + (mp_word)borrow;
            borrow = (mp_digit)((tmpw >> MP_DIGIT_BIT) & MP_MASK);
            un.dp[i + j] = (mp_digit)(tmpw & MP_MASK);
         }
      }
      /*
         D7. [Loop on j]
            Decrease j by 1;
            Test if j is not less than 0;
            If yes, go back to step D3.
       */
   }
   if (Q != NULL) {
      mp_clamp(&q);
      mp_exch(&q, Q);
      Q->sign = (neg ? MP_NEG : MP_ZPOS);
   }
   /*
      D8. [Unnormalize]
         Now (Q[m]...Q[1]Q[0]) is the desired quotient Q, and the desired remainder R may be obtained by dividing (U[n-1]...U[1]U[0]) by D.
    */
   if (R != NULL) {
      mp_clamp(&un);
      un.sign = mp_iszero(&un) ? MP_ZPOS : u->sign;
      if ((err = mp_div_2d(&un, s, R, NULL)) != MP_OKAY)        goto LTM_ERR;
   }

LTM_ERR:
   mp_clear_multi(&un, &vn, NULL);
LTM_ERR_0:
   mp_clear(&q);
   return err;
}



#else


mp_err s_mp_div_school(const mp_int *a, const mp_int *b, mp_int *c, mp_int *d)
{
   mp_int q, x, y, t1, t2;
   int n, t, i, norm;
   bool neg;
   mp_err err;

   if ((err = mp_init_size(&q, a->used + 2)) != MP_OKAY) {
      return err;
   }
   q.used = a->used + 2;

   if ((err = mp_init(&t1)) != MP_OKAY)                           goto LBL_Q;
   if ((err = mp_init(&t2)) != MP_OKAY)                           goto LBL_T1;
   if ((err = mp_init_copy(&x, a)) != MP_OKAY)                    goto LBL_T2;
   if ((err = mp_init_copy(&y, b)) != MP_OKAY)                    goto LBL_X;

   /* fix the sign */
   neg = (a->sign != b->sign);
   x.sign = y.sign = MP_ZPOS;

   /* normalize both x and y, ensure that y >= b/2, [b == 2**MP_DIGIT_BIT] */
   norm = mp_count_bits(&y) % MP_DIGIT_BIT;
   if (norm < (MP_DIGIT_BIT - 1)) {
      norm = (MP_DIGIT_BIT - 1) - norm;
      if ((err = mp_mul_2d(&x, norm, &x)) != MP_OKAY)             goto LBL_Y;
      if ((err = mp_mul_2d(&y, norm, &y)) != MP_OKAY)             goto LBL_Y;
   } else {
      norm = 0;
   }

   /* note hac does 0 based, so if used==5 then its 0,1,2,3,4, e.g. use 4 */
   n = x.used - 1;
   t = y.used - 1;

   /* while (x >= y*b**n-t) do { q[n-t] += 1; x -= y*b**{n-t} } */
   /* y = y*b**{n-t} */
   if ((err = mp_lshd(&y, n - t)) != MP_OKAY)                     goto LBL_Y;

   while (mp_cmp(&x, &y) != MP_LT) {
      ++(q.dp[n - t]);
      if ((err = mp_sub(&x, &y, &x)) != MP_OKAY)                  goto LBL_Y;
   }

   /* reset y by shifting it back down */
   mp_rshd(&y, n - t);

   /* step 3. for i from n down to (t + 1) */
   for (i = n; i >= (t + 1); i--) {
      if (i > x.used) {
         continue;
      }

      /* step 3.1 if xi == yt then set q{i-t-1} to b-1,
       * otherwise set q{i-t-1} to (xi*b + x{i-1})/yt */
      if (x.dp[i] == y.dp[t]) {
         q.dp[(i - t) - 1] = ((mp_digit)1 << (mp_digit)MP_DIGIT_BIT) - (mp_digit)1;
      } else {
         mp_word tmp;
         tmp = (mp_word)x.dp[i] << (mp_word)MP_DIGIT_BIT;
         tmp |= (mp_word)x.dp[i - 1];
         tmp /= (mp_word)y.dp[t];
         if (tmp > (mp_word)MP_MASK) {
            tmp = MP_MASK;
         }
         q.dp[(i - t) - 1] = (mp_digit)(tmp & (mp_word)MP_MASK);
      }

      /* while (q{i-t-1} * (yt * b + y{t-1})) >
               xi * b**2 + xi-1 * b + xi-2

         do q{i-t-1} -= 1;
      */
      q.dp[(i - t) - 1] = (q.dp[(i - t) - 1] + 1uL) & (mp_digit)MP_MASK;
      do {
         q.dp[(i - t) - 1] = (q.dp[(i - t) - 1] - 1uL) & (mp_digit)MP_MASK;

         /* find left hand */
         mp_zero(&t1);
         t1.dp[0] = ((t - 1) < 0) ? 0u : y.dp[t - 1];
         t1.dp[1] = y.dp[t];
         t1.used = 2;
         if ((err = mp_mul_d(&t1, q.dp[(i - t) - 1], &t1)) != MP_OKAY)   goto LBL_Y;

         /* find right hand */
         t2.dp[0] = ((i - 2) < 0) ? 0u : x.dp[i - 2];
         t2.dp[1] = x.dp[i - 1]; /* i >= 1 always holds */
         t2.dp[2] = x.dp[i];
         t2.used = 3;
      } while (mp_cmp_mag(&t1, &t2) == MP_GT);

      /* step 3.3 x = x - q{i-t-1} * y * b**{i-t-1} */
      if ((err = mp_mul_d(&y, q.dp[(i - t) - 1], &t1)) != MP_OKAY)       goto LBL_Y;
      if ((err = mp_lshd(&t1, (i - t) - 1)) != MP_OKAY)                  goto LBL_Y;
      if ((err = mp_sub(&x, &t1, &x)) != MP_OKAY)                        goto LBL_Y;

      /* if x < 0 then { x = x + y*b**{i-t-1}; q{i-t-1} -= 1; } */
      if (mp_isneg(&x)) {
         if ((err = mp_copy(&y, &t1)) != MP_OKAY)                        goto LBL_Y;
         if ((err = mp_lshd(&t1, (i - t) - 1)) != MP_OKAY)               goto LBL_Y;
         if ((err = mp_add(&x, &t1, &x)) != MP_OKAY)                     goto LBL_Y;

         q.dp[(i - t) - 1] = (q.dp[(i - t) - 1] - 1uL) & MP_MASK;
      }
   }

   /* now q is the quotient and x is the remainder
    * [which we have to normalize]
    */

   /* get sign before writing to c */
   x.sign = mp_iszero(&x) ? MP_ZPOS : a->sign;

   if (c != NULL) {
      mp_clamp(&q);
      mp_exch(&q, c);
      c->sign = (neg ? MP_NEG : MP_ZPOS);
   }

   if (d != NULL) {
      if ((err = mp_div_2d(&x, norm, &x, NULL)) != MP_OKAY)       goto LBL_Y;
      mp_exch(&x, d);
   }

LBL_Y:
   mp_clear(&y);
LBL_X:
   mp_clear(&x);
LBL_T2:
   mp_clear(&t2);
LBL_T1:
   mp_clear(&t1);
LBL_Q:
   mp_clear(&q);
   return err;
}
#endif

#endif
