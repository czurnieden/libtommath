#include "tommath_private.h"
#ifdef S_MP_SQR_TOOM_4_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



/*
   This file contains code from J. Arndt's book  "Matters Computational"
   and the accompanying FXT-library with permission of the author.
*/

/*
   Setup and interpolation from

     Chung, Jaewook, and M. Anwar Hasan. "Asymmetric squaring formulae."
     18th IEEE Symposium on Computer Arithmetic (ARITH'07). IEEE, 2007.

*/


mp_err s_mp_sqr_toom_4(const mp_int *a, mp_int *c)
{
   mp_int S0, S1, S2, S3, S4, S5, S6, S7, a0, a1, a2, a3;
   mp_err err;
   int B;

   if ((err = mp_init_multi(&S0, &S1, &S2, &S3, &S4, &S5, &S6, &S7, NULL)) != MP_OKAY) {
      return err;
   }

   B = (a->used) / 4;

   if ((err = mp_init_size(&a0, B)) != MP_OKAY)                                                          goto LTM_ERRa0;
   if ((err = mp_init_size(&a1, B)) != MP_OKAY)                                                          goto LTM_ERRa1;
   if ((err = mp_init_size(&a2, B)) != MP_OKAY)                                                          goto LTM_ERRa2;
   if ((err = mp_init_size(&a3, a->used - 3 * B)) != MP_OKAY)                                            goto LTM_ERRa3;

   /** A = a3*x^3 + a2*x^2 + a1*x + a0; */
   a0.used = a1.used = a2.used = B;
   a3.used = a->used - 3 * B;
   s_mp_copy_digs(a0.dp, a->dp, a0.used);
   s_mp_copy_digs(a1.dp, a->dp + B, a1.used);
   s_mp_copy_digs(a2.dp, a->dp + 2 * B, a2.used);
   s_mp_copy_digs(a3.dp, a->dp + 3 * B, a3.used);
   mp_clamp(&a0);
   mp_clamp(&a1);
   mp_clamp(&a2);
   mp_clamp(&a3);


#if 0
   /* This is the "naked" version and for testing purposes only! */

   /** S2 = a1 * a0; */
   if ((err = mp_mul(&a1, &a0, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2<<1; */
   if ((err = mp_mul_2d(&S2, 1, &S2)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = a1^2; */
   if ((err = mp_sqr(&a1, &S1)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S3 = a2 * a0; */
   if ((err = mp_mul(&a2, &a0, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3<<1; */
   if ((err = mp_mul_2d(&S3, 1, &S3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S3 + S1; */
   if ((err = mp_add(&S3, &S1, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = a2 * a1; */
   if ((err = mp_mul(&a2, &a1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<1; */
   if ((err = mp_mul_2d(&S1, 1, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S4 = a3 * a0; */
   if ((err = mp_mul(&a3, &a0, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4<<1; */
   if ((err = mp_mul_2d(&S4, 1, &S4)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S4 = S4 + S1; */
   if ((err = mp_add(&S4, &S1, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = a2^2; */
   if ((err = mp_sqr(&a2, &S1)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S5 = a3 * a1; */
   if ((err = mp_mul(&a3, &a1, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5<<1; */
   if ((err = mp_mul_2d(&S5, 1, &S5)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S5 = S5 + S1; */
   if ((err = mp_add(&S5, &S1, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = a3 * a2; */
   if ((err = mp_mul(&a3, &a2, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6<<1; */
   if ((err = mp_mul_2d(&S6, 1, &S6)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = a3^2; */
   if ((err = mp_sqr(&a3, &S7)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S1 = a0^2; */
   if ((err = mp_sqr(&a0, &S1)) != MP_OKAY)                                                              goto LTM_ERR;

   /** P = S1 * x^0 +  S2 * x^1 + S3 * x^2 + S4 * x^3 + S5 * x^4 + S6 * x^5 + S7 * x^6 */

   if ((err = mp_lshd(&S7, 6 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_lshd(&S6, 5 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S7, &S6, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S5, 4 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S7, &S5, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S4, 3 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S7, &S4, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S3, 2 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S7, &S3, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S2, 1 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S7, &S2, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S7, &S1, &S7)) != MP_OKAY)                                                         goto LTM_ERR;

   mp_exch(&S7, c);

   /** P - A^2*/
#endif


   /** S1 = a0 * a1; */
   if ((err = mp_mul(&a0, &a1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 * 2 ; */
   if ((err = mp_mul_2(&S1,&S1)) != MP_OKAY)                                                             goto LTM_ERR;
   /** c = a0 - a2; */
   if ((err = mp_sub(&a0, &a2, c)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S3 = a1 - a3; */
   if ((err = mp_sub(&a1, &a3, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = c + S3; */
   if ((err = mp_add(c, &S3, &S2)) != MP_OKAY)                                                           goto LTM_ERR;
   /** t = c - S3; */
   if ((err = mp_sub(c, &S3, &S0)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S2 = S2 * t; */
   if ((err = mp_mul(&S2, &S0, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = c * S3; */
   if ((err = mp_mul(c, &S3, &S0)) != MP_OKAY)                                                           goto LTM_ERR;
   /** t = t * 2 ; */
   if ((err = mp_mul_2(&S0, &S0)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S3 = a0 + a1; */
   if ((err = mp_add(&a0, &a1, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 + a2; */
   if ((err = mp_add(&S3, &a2, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 + a3; */
   if ((err = mp_add(&S3, &a3, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3^2; */
   if ((err = mp_sqr(&S3, &S3)) != MP_OKAY)                                                              goto LTM_ERR;

   /** c = a3 * a2; */
   if ((err = mp_mul(&a3, &a2, c)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = c * 2 ; */
   if ((err = mp_mul_2(c, c)) != MP_OKAY)                                                                goto LTM_ERR;
   /** a1 = S2 + S3; */
   if ((err = mp_add(&S2, &S3, &a1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a1 = a1 + t; */
   if ((err = mp_add(&a1, &S0, &a1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a1 = a1 / 2 ; */
   if ((err = mp_div_2(&a1, &a1)) != MP_OKAY)                                                            goto LTM_ERR;

   /** a2 = S1 + c; */
   if ((err = mp_add(&S1, c, &a2)) != MP_OKAY)                                                           goto LTM_ERR;
   /** a1 = a1 - a2; */
   if ((err = mp_sub(&a1, &a2, &a1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = a2 - t; */
   if ((err = mp_sub(&a2, &S0, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = a1 - S2; */
   if ((err = mp_sub(&a1, &S2, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a0 = a0^2; */
   if ((err = mp_sqr(&a0, &a0)) != MP_OKAY)                                                              goto LTM_ERR;
   /** t = a1 - a0; */
   if ((err = mp_sub(&a1, &a0, &S0)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a3 = a3^2; */
   if ((err = mp_sqr(&a3, &a3)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S2 = S2 - a3; */
   if ((err = mp_sub(&S2, &a3, &S2)) != MP_OKAY)                                                         goto LTM_ERR;

   /** P = a3 * x^6 + c * x^5 + t * x^4 + S3 * x^3 + S2 * x^2 + S1 * x + a0; */
   if ((err = mp_lshd(&a3, 6 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_lshd(c, 5 * B)) != MP_OKAY)                                                             goto LTM_ERR;
   if ((err = mp_add(&a3, c, &a3)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_lshd(&S0, 4 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&a3, &S0, &S0)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S3, 3 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S3, &S0, &S0)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S2, 2 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S2, &S0, &S0)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S1, 1 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S0, &S0)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&a0, &S0, c)) != MP_OKAY)                                                           goto LTM_ERR;

   /** A^2 - P */

LTM_ERR:
   mp_clear(&a3);
LTM_ERRa3:
   mp_clear(&a2);
LTM_ERRa2:
   mp_clear(&a1);
LTM_ERRa1:
   mp_clear(&a0);
LTM_ERRa0:
   mp_clear_multi(&S0, &S1, &S2, &S3, &S4, &S5, &S6, &S7, NULL);
   return err;
}


#endif
