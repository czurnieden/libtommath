#include "tommath_private.h"
#ifdef S_MP_MUL_TOOM_4_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



/*
   This file contains code from J. Arndt's book  "Matters Computational"
   and the accompanying FXT-library with permission of the author.
*/

/*
    Bodrato, Marco, and Alberto Zanoni. "Integer and polynomial multiplication:
    Towards optimal Toom-Cook matrices." Proceedings of the 2007 international
    symposium on Symbolic and algebraic computation. ACM, 2007.
*/
mp_err s_mp_mul_toom_4(const mp_int *a, const mp_int *b, mp_int *c)
{
   mp_int S1, S2, S3, S4, S5, S6, S7;
   mp_int a0, a1, a2, a3;
   mp_int b0, b1, b2, b3;

   mp_err err = MP_OKAY;
   int B;

   B = MP_MIN(a->used, b->used) / 4;

   if ((err = mp_init_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, NULL)) != MP_OKAY) {
      return err;
   }

   if ((err = mp_init_size(&a0, B)) != MP_OKAY)                                                          goto LTM_ERRa0;
   if ((err = mp_init_size(&a1, B)) != MP_OKAY)                                                          goto LTM_ERRa1;
   if ((err = mp_init_size(&a2, B)) != MP_OKAY)                                                          goto LTM_ERRa2;
   if ((err = mp_init_size(&a3, a->used - 3 * B)) != MP_OKAY)                                            goto LTM_ERRa3;

   /** a = a3*B^3 + a2*B^2 + a1*B + a0; */
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


   if ((err = mp_init_size(&b0, B)) != MP_OKAY)                                                          goto LTM_ERRb0;
   if ((err = mp_init_size(&b1, B)) != MP_OKAY)                                                          goto LTM_ERRb1;
   if ((err = mp_init_size(&b2, B)) != MP_OKAY)                                                          goto LTM_ERRb2;
   if ((err = mp_init_size(&b3, b->used - 3 * B)) != MP_OKAY)                                            goto LTM_ERRb3;

   /** b = b3*B^3 + b2*B^2 + b1*B + b0; */
   b0.used = b1.used = b2.used = B;
   b3.used = b->used - 3 * B;
   s_mp_copy_digs(b0.dp, b->dp, b0.used);
   s_mp_copy_digs(b1.dp, b->dp + B, b1.used);
   s_mp_copy_digs(b2.dp, b->dp + 2 * B, b2.used);
   s_mp_copy_digs(b3.dp, b->dp + 3 * B, b3.used);
   mp_clamp(&b0);
   mp_clamp(&b1);
   mp_clamp(&b2);
   mp_clamp(&b3);


   /** \\S2 = (8*a3 + 4*a2 + 2*a1 + a0)*(8*b3 + 4*b2 + 2*b1 + b0); */
   /** S1 = a3 * 2; */
   if ((err = mp_mul_2(&a3, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = S1 + a2; */
   if ((err = mp_add(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 * 2; */
   if ((err = mp_mul_2(&S1, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = S1 + a1; */
   if ((err = mp_add(&S1, &a1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 * 2; */
   if ((err = mp_mul_2(&S1, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = S1 + a0; */
   if ((err = mp_add(&S1, &a0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;

   /** S2 = b3 * 2; */
   if ((err = mp_mul_2(&b3, &S2)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S2 = S2 + b2; */
   if ((err = mp_add(&S2, &b2, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 * 2; */
   if ((err = mp_mul_2(&S2, &S2)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S2 = S2 + b1; */
   if ((err = mp_add(&S2, &b1, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 * 2; */
   if ((err = mp_mul_2(&S2, &S2)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S2 = S2 + b0; */
   if ((err = mp_add(&S2, &b0, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S1 * S2; */
   if ((err = mp_mul(&S1, &S2, &S2)) != MP_OKAY)                                                         goto LTM_ERR;

   /** \\S3 = (+a3 + a2 + a1 + a0)*(+b3 + b2 + b1 + b0); */
   /** S1 = a3 + a2; */
   if ((err = mp_add(&a3, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + a1; */
   if ((err = mp_add(&S1, &a1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + a0; */
   if ((err = mp_add(&S1, &a0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = b3 + b2; */
   if ((err = mp_add(&b3, &b2, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 + b1; */
   if ((err = mp_add(&S3, &b1, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 + b0; */
   if ((err = mp_add(&S3, &b0, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S1 * S3; */
   if ((err = mp_mul(&S1, &S3, &S3)) != MP_OKAY)                                                         goto LTM_ERR;

   /** \\S4 = (-a3 + a2 - a1 + a0) * ( -b3 + b2 - b1 + b0); */
   /** \\S4 = (a2 + a0 -a1 - a3) * (b2 + b0 -b1 - b3) */
   /** S1 = a2 + a0; */
   if ((err = mp_add(&a2, &a0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 - a1; */
   if ((err = mp_sub(&S1, &a1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 - a3; */
   if ((err = mp_sub(&S1, &a3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = b2 + b0; */
   if ((err = mp_add(&b2, &b0, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - b1; */
   if ((err = mp_sub(&S4, &b1, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - b3; */
   if ((err = mp_sub(&S4, &b3, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S1 * S4; */
   if ((err = mp_mul(&S1, &S4, &S4)) != MP_OKAY)                                                         goto LTM_ERR;

   /** \\S5 = (+8*a0+4*a1+2*a2+a3)*(+8*b0+4*b1+2*b2+b3); */
   /** S1 = a0 * 2; */
   if ((err = mp_mul_2(&a0, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = S1 + a1; */
   if ((err = mp_add(&S1, &a1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 * 2; */
   if ((err = mp_mul_2(&S1, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = S1 + a2; */
   if ((err = mp_add(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 * 2; */
   if ((err = mp_mul_2(&S1, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = S1 + a3; */
   if ((err = mp_add(&S1, &a3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;

   /** S5 = b0 * 2; */
   if ((err = mp_mul_2(&b0, &S5)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = S5 + b1; */
   if ((err = mp_add(&S5, &b1, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 * 2; */
   if ((err = mp_mul_2(&S5, &S5)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = S5 + b2; */
   if ((err = mp_add(&S5, &b2, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 * 2; */
   if ((err = mp_mul_2(&S5, &S5)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = S5 + b3; */
   if ((err = mp_add(&S5, &b3, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S1 * S5; */
   if ((err = mp_mul(&S1, &S5, &S5)) != MP_OKAY)                                                         goto LTM_ERR;


   /** \\S6 = (-8*a0+4*a1-2*a2+a3)*(-8*b0+4*b1-2*b2+b3); */
   /** S1 = a0 * 2; */
   if ((err = mp_mul_2(&a0, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = -S1; */
   S1.sign = (S1.sign == MP_ZPOS)? MP_NEG: MP_ZPOS;
   /** S1 = S1 + a1; */
   if ((err = mp_add(&S1, &a1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 * 2; */
   if ((err = mp_mul_2(&S1, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = S1 - a2; */
   if ((err = mp_sub(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 * 2; */
   if ((err = mp_mul_2(&S1, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = S1 + a3; */
   if ((err = mp_add(&S1, &a3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;

   /** S6 = b0 * 2; */
   if ((err = mp_mul_2(&b0, &S6)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S6 = -S6; */
   S6.sign = (S6.sign == MP_ZPOS)? MP_NEG: MP_ZPOS;
   /** S6 = S6 + b1; */
   if ((err = mp_add(&S6, &b1, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 * 2; */
   if ((err = mp_mul_2(&S6, &S6)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S6 = S6 - b2; */
   if ((err = mp_sub(&S6, &b2, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 * 2; */
   if ((err = mp_mul_2(&S6, &S6)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S6 = S6 + b3; */
   if ((err = mp_add(&S6, &b3, &S6)) != MP_OKAY)                                                         goto LTM_ERR;

   /** S6 = S1 * S6; */
   if ((err = mp_mul(&S1, &S6, &S6)) != MP_OKAY)                                                         goto LTM_ERR;

   /** S1 = a3 * b3; */
   if ((err = mp_mul(&a3, &b3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;

   /** S7 = a0 * b0; */
   if ((err = mp_mul(&a0, &b0, &S7)) != MP_OKAY)                                                         goto LTM_ERR;

   /** S2 = S2 + S5; */
   if ((err = mp_add(&S2, &S5, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - S3; */
   if ((err = mp_sub(&S4, &S3, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - S5; */
   if ((err = mp_sub(&S6, &S5, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 /2; */
   if ((err = mp_div_2d(&S4, 1, &S4, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S5 = S5 - S1; */
   if ((err = mp_sub(&S5, &S1, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 - (64*S7); */
   if ((err = mp_mul_2d(&S7, 6, c)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_sub(&S5, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S3 = S3 + S4; */
   if ((err = mp_add(&S3, &S4, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 *2 ; */
   if ((err = mp_mul_2d(&S5, 1, &S5)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S5 = S5 + S6; */
   if ((err = mp_add(&S5, &S6, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 - (65*S3) ; */
   if ((err = mp_mul_d(&S3, 65, c)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_sub(&S2, c, &S2)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S3 = S3 - S1; */
   if ((err = mp_sub(&S3, &S1, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 - S7; */
   if ((err = mp_sub(&S3, &S7, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = -S4; */
   if ((err = mp_neg(&S4, &S4)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S6 = -S6; */
   if ((err = mp_neg(&S6, &S6)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S2 = S2 + (45*S3) ; */
   if ((err = mp_mul_d(&S3, 45, c)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_add(&S2, c, &S2)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S5 = S5 - (8*S3) ; */
   if ((err = mp_mul_2d(&S3, 3, c)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_sub(&S5, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S5 = S5/24; */
   if ((err = mp_div_d(&S5, 24, &S5, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S6 = S6 - S2; */
   if ((err = mp_sub(&S6, &S2, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 - (16*S4) ; */
   if ((err = mp_mul_2d(&S4, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_sub(&S2, c, &S2)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S2 = S2/18; */
   if ((err = mp_div_d(&S2, 18, &S2, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S3 = S3 - S5; */
   if ((err = mp_sub(&S3, &S5, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - S2; */
   if ((err = mp_sub(&S4, &S2, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 + (30*S2) ; */
   if ((err = mp_mul_d(&S2, 30, c)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_add(&S6, c, &S6)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S6 = S6/60; */
   if ((err = mp_div_d(&S6, 60, &S6, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S2 = S2 - S6; */
   if ((err = mp_sub(&S2, &S6, &S2)) != MP_OKAY)                                                         goto LTM_ERR;

   /** S = S1*B^6 + S2*B^5 + S3*B^4 + S4*B^3 + S5*B^2 + S6*B + S7; */
   if ((err = mp_lshd(&S1, 6 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_lshd(&S2, 5 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S3, 4 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S4, 3 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S5, 2 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S5, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S6, 1 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S1, &S7, c)) != MP_OKAY)                                                           goto LTM_ERR;

   /** S - a*b */

LTM_ERR:
   mp_clear(&b3);
LTM_ERRb3:
   mp_clear(&b2);
LTM_ERRb2:
   mp_clear(&b1);
LTM_ERRb1:
   mp_clear(&b0);
LTM_ERRb0:
   mp_clear(&a3);
LTM_ERRa3:
   mp_clear(&a2);
LTM_ERRa2:
   mp_clear(&a1);
LTM_ERRa1:
   mp_clear(&a0);
LTM_ERRa0:
   mp_clear_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, NULL);
   return err;
}



#endif
