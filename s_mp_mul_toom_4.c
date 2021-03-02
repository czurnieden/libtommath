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
   mp_int S1, S2, S3, S4, S5, S6, S7, t;
   mp_int a0, a1, a2, a3;
   mp_int b0, b1, b2, b3;

   mp_err err = MP_OKAY;
   int B;

   B = MP_MIN(a->used, b->used) / 4;

   if ((err = mp_init_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, &t, NULL)) != MP_OKAY) {
      return err;
   }

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


   if ((err = mp_init_size(&b0, B)) != MP_OKAY)                                                          goto LTM_ERRb0;
   if ((err = mp_init_size(&b1, B)) != MP_OKAY)                                                          goto LTM_ERRb1;
   if ((err = mp_init_size(&b2, B)) != MP_OKAY)                                                          goto LTM_ERRb2;
   if ((err = mp_init_size(&b3, b->used - 3 * B)) != MP_OKAY)                                            goto LTM_ERRb3;

   /** B = b3*x^3 + b2*x^2 + b1*x + b0; */
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


#if 0
   /** S5 = a3 * 2; */
   if ((err = mp_mul_2(&a3, &S5)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = S5 + a2; */
   if ((err = mp_add(&S5, &a2, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 * 2; */
   if ((err = mp_mul_2(&S5, &S5)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = S5 + a1; */
   if ((err = mp_add(&S5, &a1, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 * 2; */
   if ((err = mp_mul_2(&S5, &S5)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = S5 + a0; */
   if ((err = mp_add(&S5, &a0, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = b3 * 2; */
   if ((err = mp_mul_2(&b3, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = S1 + b2; */
   if ((err = mp_add(&S1, &b2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 * 2; */
   if ((err = mp_mul_2(&S1, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = S1 + b1; */
   if ((err = mp_add(&S1, &b1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 * 2; */
   if ((err = mp_mul_2(&S1, &S1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S1 = S1 + b0; */
   if ((err = mp_add(&S1, &b0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S5 * S1; */
   if ((err = mp_mul(&S5, &S1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = a3 + a2; */
   if ((err = mp_add(&a3, &a2, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 + a1; */
   if ((err = mp_add(&S5, &a1, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 + a0; */
   if ((err = mp_add(&S5, &a0, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = b3 + b2; */
   if ((err = mp_add(&b3, &b2, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 + b1; */
   if ((err = mp_add(&S2, &b1, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 + b0; */
   if ((err = mp_add(&S2, &b0, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S5 * S2; */
   if ((err = mp_mul(&S5, &S2, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = a2 + a0; */
   if ((err = mp_add(&a2, &a0, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 - a1; */
   if ((err = mp_sub(&S5, &a1, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 - a3; */
   if ((err = mp_sub(&S5, &a3, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = b2 + b0; */
   if ((err = mp_add(&b2, &b0, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 - b1; */
   if ((err = mp_sub(&S3, &b1, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 - b3; */
   if ((err = mp_sub(&S3, &b3, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S5 * S3; */
   if ((err = mp_mul(&S5, &S3, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** c = a0 * 2; */
   if ((err = mp_mul_2(&a0, c)) != MP_OKAY)                                                              goto LTM_ERR;
   /** c = c + a1; */
   if ((err = mp_add(c, &a1, c)) != MP_OKAY)                                                             goto LTM_ERR;
   /** c = c * 2; */
   if ((err = mp_mul_2(c, c)) != MP_OKAY)                                                                goto LTM_ERR;
   /** c = c + a2; */
   if ((err = mp_add(c, &a2, c)) != MP_OKAY)                                                             goto LTM_ERR;
   /** c = c * 2; */
   if ((err = mp_mul_2(c, c)) != MP_OKAY)                                                                goto LTM_ERR;
   /** c = c + a3; */
   if ((err = mp_add(c, &a3, c)) != MP_OKAY)                                                             goto LTM_ERR;
   /** S5 = a0 * 2; */
   if ((err = mp_mul_2(&a0, &S5)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = a1 - S5; */
   if ((err = mp_sub(&a1, &S5, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a1 = b0 * 2; */
   if ((err = mp_mul_2(&b0, &a1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** a1 = a1 + b1; */
   if ((err = mp_add(&a1, &b1, &a1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a1 = a1 * 2; */
   if ((err = mp_mul_2(&a1, &a1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** a1 = a1 + b2; */
   if ((err = mp_add(&a1, &b2, &a1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a1 = a1 * 2; */
   if ((err = mp_mul_2(&a1, &a1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** a1 = a1 + b3; */
   if ((err = mp_add(&a1, &b3, &a1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a1 = c * a1; */
   if ((err = mp_mul(c, &a1, &a1)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S5 = S5 * 2; */
   if ((err = mp_mul_2(&S5, &S5)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = S5 - a2; */
   if ((err = mp_sub(&S5, &a2, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 * 2; */
   if ((err = mp_mul_2(&S5, &S5)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = S5 + a3; */
   if ((err = mp_add(&S5, &a3, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a2 = b0 * 2; */
   if ((err = mp_mul_2(&b0, &a2)) != MP_OKAY)                                                            goto LTM_ERR;
   /** a2 = -a2; */
   a2.sign = ((!mp_iszero(&a2) && !mp_isneg(&a2)) ? MP_NEG : MP_ZPOS);
   /** a2 = a2 + b1; */
   if ((err = mp_add(&a2, &b1, &a2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a2 = a2 * 2; */
   if ((err = mp_mul_2(&a2, &a2)) != MP_OKAY)                                                            goto LTM_ERR;
   /** a2 = a2 - b2; */
   if ((err = mp_sub(&a2, &b2, &a2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a2 = a2 * 2; */
   if ((err = mp_mul_2(&a2, &a2)) != MP_OKAY)                                                            goto LTM_ERR;
   /** a2 = a2 + b3; */
   if ((err = mp_add(&a2, &b3, &a2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a2 = S5 * a2; */
   if ((err = mp_mul(&S5, &a2, &a2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = a3 * b3; */
   if ((err = mp_mul(&a3, &b3, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** b3 = a0 * b0; */
   if ((err = mp_mul(&a0, &b0, &b3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + a1; */
   if ((err = mp_add(&S1, &a1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 - S2; */
   if ((err = mp_sub(&S3, &S2, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a2 = a2 - a1; */
   if ((err = mp_sub(&a2, &a1, &a2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 / 2; */
   if ((err = mp_div_2(&S3, &S3)) != MP_OKAY)                                                            goto LTM_ERR;
   /** a1 = a1 - S5; */
   if ((err = mp_sub(&a1, &S5, &a1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a3 = b3<<6; */
   if ((err = mp_mul_2d(&b3, 6, &a3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** a1 = a1 - a3; */
   if ((err = mp_sub(&a1, &a3, &a1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 + S3; */
   if ((err = mp_add(&S2, &S3, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a1 = a1 * 2; */
   if ((err = mp_mul_2(&a1, &a1)) != MP_OKAY)                                                            goto LTM_ERR;
   /** a1 = a1 + a2; */
   if ((err = mp_add(&a1, &a2, &a1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a3 = S2 * 65; */
   if ((err = mp_mul_d(&S2, 65, &a3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 - a3; */
   if ((err = mp_sub(&S1, &a3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 - S5; */
   if ((err = mp_sub(&S2, &S5, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 - b3; */
   if ((err = mp_sub(&S2, &b3, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = -S3; */
   S3.sign = ((!mp_iszero(&S3) && !mp_isneg(&S3)) ? MP_NEG : MP_ZPOS);
   /** a2 = -a2; */
   a2.sign = ((!mp_iszero(&a2) && !mp_isneg(&a2)) ? MP_NEG : MP_ZPOS);
   /** a3 = S2 * 45; */
   if ((err = mp_mul_d(&S2, 45, &a3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a3; */
   if ((err = mp_add(&S1, &a3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a3 = S2<<3; */
   if ((err = mp_mul_2d(&S2, 3, &a3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** a1 = a1 - a3; */
   if ((err = mp_sub(&a1, &a3, &a1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a1 = a1 / 24; */
   if ((err = mp_div_d(&a1, 24, &a1, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** a2 = a2 - S1; */
   if ((err = mp_sub(&a2, &S1, &a2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a3 = S3<<4; */
   if ((err = mp_mul_2d(&S3, 4, &a3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 - a3; */
   if ((err = mp_sub(&S1, &a3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 / 18; */
   if ((err = mp_div_d(&S1, 18, &S1, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S2 = S2 - a1; */
   if ((err = mp_sub(&S2, &a1, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 - S1; */
   if ((err = mp_sub(&S3, &S1, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a3 = S1 * 30; */
   if ((err = mp_mul_d(&S1, 30, &a3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** a2 = a2 + a3; */
   if ((err = mp_add(&a2, &a3, &a2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** a2 = a2 / 60; */
   if ((err = mp_div_d(&a2, 60, &a2, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S1 = S1 - a2; */
   if ((err = mp_sub(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;

   /** P = S5 * x^6 + S1 * x^5 + S2 * x^4 + S3 * x^3 + a1 * x^2 + a2 * x + b3;*/
   if ((err = mp_lshd(&S5, 6 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_lshd(&S1, 5 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S5, &S1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S2, 4 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S3, 3 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&a1, 2 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &a1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&a2, 1 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S1, &b3, c)) != MP_OKAY)                                                           goto LTM_ERR;

   /** S - a*b */
#endif

   /** S1 = a2<<2; */
   if ((err = mp_mul_2d(&a2, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a0; */
   if ((err = mp_add(&S1, &a0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = a3<<2; */
   if ((err = mp_mul_2d(&a3, 2, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = S7 + a1; */
   if ((err = mp_add(&S7, &a1, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7<<1; */
   if ((err = mp_mul_2d(&S7, 1, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S1 - S7; */
   if ((err = mp_sub(&S1, &S7, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S1 + S7; */
   if ((err = mp_add(&S1, &S7, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = b2<<2; */
   if ((err = mp_mul_2d(&b2, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b0; */
   if ((err = mp_add(&S1, &b0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = b3<<2; */
   if ((err = mp_mul_2d(&b3, 2, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = S7 + b1; */
   if ((err = mp_add(&S7, &b1, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7<<1; */
   if ((err = mp_mul_2d(&S7, 1, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S4 = S1 - S7; */
   if ((err = mp_sub(&S1, &S7, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S1 + S7; */
   if ((err = mp_add(&S1, &S7, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S6 * S5; */
   if ((err = mp_mul(&S6, &S5, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 * S4; */
   if ((err = mp_mul(&S3, &S4, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = a0 + a2; */
   if ((err = mp_add(&a0, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = a1 + a3; */
   if ((err = mp_add(&a1, &a3, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S1 - S7; */
   if ((err = mp_sub(&S1, &S7, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S1 + S7; */
   if ((err = mp_add(&S1, &S7, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = b0 + b2; */
   if ((err = mp_add(&b0, &b2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = b1 + b3; */
   if ((err = mp_add(&b1, &b3, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S1 - S7; */
   if ((err = mp_sub(&S1, &S7, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S1 + S7; */
   if ((err = mp_add(&S1, &S7, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 * S7; */
   if ((err = mp_mul(&S4, &S7, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 * S6; */
   if ((err = mp_mul(&S5, &S6, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = a0<<2; */
   if ((err = mp_mul_2d(&a0, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a2; */
   if ((err = mp_add(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<1; */
   if ((err = mp_mul_2d(&S1, 1, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = a1<<2; */
   if ((err = mp_mul_2d(&a1, 2, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = S7 + a3; */
   if ((err = mp_add(&S7, &a3, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S1 - S7; */
   if ((err = mp_sub(&S1, &S7, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = b0<<2; */
   if ((err = mp_mul_2d(&b0, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b2; */
   if ((err = mp_add(&S1, &b2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<1; */
   if ((err = mp_mul_2d(&S1, 1, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = b1<<2; */
   if ((err = mp_mul_2d(&b1, 2, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = S7 + b3; */
   if ((err = mp_add(&S7, &b3, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S1 - S7; */
   if ((err = mp_sub(&S1, &S7, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 * S7; */
   if ((err = mp_mul(&S6, &S7, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = a0 * b0; */
   if ((err = mp_mul(&a0, &b0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = a3 * b3; */
   if ((err = mp_mul(&a3, &b3, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 - S2; */
   if ((err = mp_sub(&S3, &S2, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = -S3; */
   S3.sign = (S3.sign == MP_NEG)?MP_ZPOS:MP_NEG;
   /** S3 = S3>>1; */
   if ((err = mp_div_2d(&S3, 1, &S3, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S2 = S2 - S3; */
   if ((err = mp_sub(&S2, &S3, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 - S1; */
   if ((err = mp_sub(&S2, &S1, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 + S4; */
   if ((err = mp_add(&S5, &S4, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5>>1; */
   if ((err = mp_div_2d(&S5, 1, &S5, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S4 = S4 - S5; */
   if ((err = mp_sub(&S4, &S5, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - S7; */
   if ((err = mp_sub(&S6, &S7, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 + S2; */
   if ((err = mp_add(&S6, &S2, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - S3; */
   if ((err = mp_sub(&S6, &S3, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S4<<4; */
   if ((err = mp_mul_2d(&S4, 4, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 + t; */
   if ((err = mp_add(&S6, &t, &S6)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S7<<6; */
   if ((err = mp_mul_2d(&S7, 6, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S6 = S6 - t; */
   if ((err = mp_sub(&S6, &t, &S6)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S1<<6; */
   if ((err = mp_mul_2d(&S1, 6, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - t; */
   if ((err = mp_sub(&S6, &t, &S6)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S5 = S5 - S1; */
   if ((err = mp_sub(&S5, &S1, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 - S7; */
   if ((err = mp_sub(&S5, &S7, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S4<<3; */
   if ((err = mp_mul_2d(&S4, 3, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 - t; */
   if ((err = mp_sub(&S3, &t, &S3)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S5 * 20; */
   if ((err = mp_mul_d(&S5, 20, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - t; */
   if ((err = mp_sub(&S6, &t, &S6)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S6 = -S6; */
   S6.sign = (S6.sign == MP_NEG)?MP_ZPOS:MP_NEG;
   /** S6 = S6 / 18; */
   if ((err = mp_div_d(&S6, 18, &S6, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S4 = S4 - S6; */
   if ((err = mp_sub(&S4, &S6, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S6 * 6; */
   if ((err = mp_mul_d(&S6, 6, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = S3 + t; */
   if ((err = mp_add(&S3, &t, &S3)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S5<<2; */
   if ((err = mp_mul_2d(&S5, 2, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S2 = S2 / 12; */
   if ((err = mp_div_d(&S2, 12, &S2, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S3 = S3 / 30; */
   if ((err = mp_div_d(&S3, 30, &S3, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S6 = S6 - S3; */
   if ((err = mp_sub(&S6, &S3, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 - S2; */
   if ((err = mp_sub(&S5, &S2, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S2;S2 = S6;S6 = t;*/
   mp_exch(&S2, &S6);
   /** t = S3;S3 = S5;S5 = t;*/
   mp_exch(&S3, &S5);
   /** t = S5;S5 = S6;S6 = t;*/
   mp_exch(&S5, &S6);
   /** P = S1*x^0 + S2*x^1 + S3*x^2 + S4*x^3 + S5*x^4 + S6*x^5 + S7*x^6; */
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
   /** P - A*B */



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
   mp_clear_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7,  &t, NULL);
   return err;
}



#endif
