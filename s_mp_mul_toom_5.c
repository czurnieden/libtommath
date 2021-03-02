#include "tommath_private.h"
#ifdef S_MP_MUL_TOOM_5_C
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

mp_err s_mp_mul_toom_5(const mp_int *a, const mp_int *b, mp_int *c)
{
   mp_int S1, S2, S3, S4, S5, S6, S7, S8, S9;
   mp_int a0, a1, a2, a3, a4;
   mp_int b0, b1, b2, b3, b4;
   mp_err err = MP_OKAY;
   int B, sign;

   B = MP_MIN(a->used, b->used) / 5;
   if ((err =
           mp_init_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, &S8, &S9, NULL)) != MP_OKAY) {
      return err;
   }

   if ((err = mp_init_size(&a0, B)) != MP_OKAY)                                                          goto LTM_ERRa0;
   if ((err = mp_init_size(&a1, B)) != MP_OKAY)                                                          goto LTM_ERRa1;
   if ((err = mp_init_size(&a2, B)) != MP_OKAY)                                                          goto LTM_ERRa2;
   if ((err = mp_init_size(&a3, B)) != MP_OKAY)                                                          goto LTM_ERRa3;
   if ((err = mp_init_size(&a4, a->used - 4 * B)) != MP_OKAY)                                            goto LTM_ERRa4;

   /** A = a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0; */
   a0.used = a1.used = a2.used = a3.used = B;
   a4.used = a->used - 4 * B;
   s_mp_copy_digs(a0.dp, a->dp, a0.used);
   s_mp_copy_digs(a1.dp, a->dp + B, a1.used);
   s_mp_copy_digs(a2.dp, a->dp + 2 * B, a2.used);
   s_mp_copy_digs(a3.dp, a->dp + 3 * B, a3.used);
   s_mp_copy_digs(a4.dp, a->dp + 4 * B, a4.used);
   mp_clamp(&a0);
   mp_clamp(&a1);
   mp_clamp(&a2);
   mp_clamp(&a3);
   mp_clamp(&a4);

   if ((err = mp_init_size(&b0, B)) != MP_OKAY)                                                          goto LTM_ERRb0;
   if ((err = mp_init_size(&b1, B)) != MP_OKAY)                                                          goto LTM_ERRb1;
   if ((err = mp_init_size(&b2, B)) != MP_OKAY)                                                          goto LTM_ERRb2;
   if ((err = mp_init_size(&b3, B)) != MP_OKAY)                                                          goto LTM_ERRb3;
   if ((err = mp_init_size(&b4, b->used - 4 * B)) != MP_OKAY)                                            goto LTM_ERRb4;

   /** B = b4*x^4 + b3*x^3 + b2*x^2 + b1*x + b0; */
   b0.used = b1.used = b2.used = b3.used = B;
   b4.used = b->used - 4 * B;
   s_mp_copy_digs(b0.dp, b->dp, b0.used);
   s_mp_copy_digs(b1.dp, b->dp + B, b1.used);
   s_mp_copy_digs(b2.dp, b->dp + 2 * B, b2.used);
   s_mp_copy_digs(b3.dp, b->dp + 3 * B, b3.used);
   s_mp_copy_digs(b4.dp, b->dp + 4 * B, b4.used);
   mp_clamp(&b0);
   mp_clamp(&b1);
   mp_clamp(&b2);
   mp_clamp(&b3);
   mp_clamp(&b4);

   /** S1 = a4*b4; */
   if ((err = mp_mul(&a4, &b4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S2 = (a0-2*a1+4*a2-8*a3+16*a4)*(b0-2*b1+4*b2-8*b3+16*b4); */
   /** c = a1 << 1; */
   if ((err = mp_mul_2(&a1, c)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S2 = a0 - c; */
   if ((err = mp_sub(&a0, c, &S2)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a2 << 2; */
   if ((err = mp_mul_2d(&a2, 2, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S2 = S2 + c; */
   if ((err = mp_add(&S2, c, &S2)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a3 << 3; */
   if ((err = mp_mul_2d(&a3, 3, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S2 = S2 - c; */
   if ((err = mp_sub(&S2, c, &S2)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a4 << 4; */
   if ((err = mp_mul_2d(&a4, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S2 = S2 + c; */
   if ((err = mp_add(&S2, c, &S2)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b1 << 1; */
   if ((err = mp_mul_2(&b1, c)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S5 = b0 - c; */
   if ((err = mp_sub(&b0, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b2 << 2; */
   if ((err = mp_mul_2d(&b2, 2, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S5 = S5 + c; */
   if ((err = mp_add(&S5, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b3 << 3; */
   if ((err = mp_mul_2d(&b3, 3, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S5 = S5 - c; */
   if ((err = mp_sub(&S5, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b4 << 4; */
   if ((err = mp_mul_2d(&b4, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S5 = S5 + c; */
   if ((err = mp_add(&S5, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S2 = S2 * S5; */
   if ((err = mp_mul(&S2, &S5, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S5 = (a0+2*a1+4*a2+8*a3+16*a4)*(b0+2*b1+4*b2+8*b3+16*b4); */
   /** c = a1 << 1; */
   if ((err = mp_mul_2(&a1, c)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S5 = a0 + c; */
   if ((err = mp_add(&a0, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a2 << 2; */
   if ((err = mp_mul_2d(&a2, 2, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S5 = S5 + c; */
   if ((err = mp_add(&S5, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a3 << 3; */
   if ((err = mp_mul_2d(&a3, 3, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S5 = S5 + c; */
   if ((err = mp_add(&S5, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a4 << 4; */
   if ((err = mp_mul_2d(&a4, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S5 = S5 + c; */
   if ((err = mp_add(&S5, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b1 << 1; */
   if ((err = mp_mul_2(&b1, c)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S3 = b0 + c; */
   if ((err = mp_add(&b0, c, &S3)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b2 << 2; */
   if ((err = mp_mul_2d(&b2, 2, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = S3 + c; */
   if ((err = mp_add(&S3, c, &S3)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b3 << 3; */
   if ((err = mp_mul_2d(&b3, 3, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = S3 + c; */
   if ((err = mp_add(&S3, c, &S3)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b4 << 4; */
   if ((err = mp_mul_2d(&b4, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = S3 + c; */
   if ((err = mp_add(&S3, c, &S3)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S5 = S5 * S3; */
   if ((err = mp_mul(&S5, &S3, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S3 = (a4+2*a3+4*a2+8*a1+16*a0)*(b4+2*b3+4*b2+8*b1+16*b0); */

   /** c = a3 << 1; */
   if ((err = mp_mul_2(&a3, c)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S3 = a4 + c; */
   if ((err = mp_add(&a4, c, &S3)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a2 << 2; */
   if ((err = mp_mul_2d(&a2, 2, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = S3 + c; */
   if ((err = mp_add(&S3, c, &S3)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a1 << 3; */
   if ((err = mp_mul_2d(&a1, 3, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = S3 + c; */
   if ((err = mp_add(&S3, c, &S3)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a0 << 4; */
   if ((err = mp_mul_2d(&a0, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = S3 + c; */
   if ((err = mp_add(&S3, c, &S3)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b3 << 1; */
   if ((err = mp_mul_2(&b3, c)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S8 = b4 + c; */
   if ((err = mp_add(&b4, c, &S8)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b2 << 2; */
   if ((err = mp_mul_2d(&b2, 2, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S8 = S8 + c; */
   if ((err = mp_add(&S8, c, &S8)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b1 << 3; */
   if ((err = mp_mul_2d(&b1, 3, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S8 = S8 + c; */
   if ((err = mp_add(&S8, c, &S8)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b0 << 4; */
   if ((err = mp_mul_2d(&b0, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S8 = S8 + c; */
   if ((err = mp_add(&S8, c, &S8)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S3 = S3 * S8; */
   if ((err = mp_mul(&S3, &S8, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = (a4-2*a3+4*a2-8*a1+16*a0)*(b4-2*b3+4*b2-8*b1+16*b0); */
   /** c = a3 << 1; */
   if ((err = mp_mul_2(&a3, c)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S8 = a4 - c; */
   if ((err = mp_sub(&a4, c, &S8)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a2 << 2; */
   if ((err = mp_mul_2d(&a2, 2, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S8 = S8 + c; */
   if ((err = mp_add(&S8, c, &S8)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a1 << 3; */
   if ((err = mp_mul_2d(&a1, 3, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S8 = S8 - c; */
   if ((err = mp_sub(&S8, c, &S8)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a0 << 4; */
   if ((err = mp_mul_2d(&a0, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S8 = S8 + c; */
   if ((err = mp_add(&S8, c, &S8)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b3 << 1; */
   if ((err = mp_mul_2(&b3, c)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S4 = b4 - c; */
   if ((err = mp_sub(&b4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b2 << 2; */
   if ((err = mp_mul_2d(&b2, 2, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 + c; */
   if ((err = mp_add(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b1 << 3; */
   if ((err = mp_mul_2d(&b1, 3, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 - c; */
   if ((err = mp_sub(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b0 << 4; */
   if ((err = mp_mul_2d(&b0, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 + c; */
   if ((err = mp_add(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S8 = S8 * S4; */
   if ((err = mp_mul(&S8, &S4, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S4 = (a0+4*a1+16*a2+64*a3+256*a4)*(b0+4*b1+16*b2+64*b3+256*b4); */
   /** c = a1 << 2; */
   if ((err = mp_mul_2d(&a1, 2, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = a0 + c; */
   if ((err = mp_add(&a0, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a2 << 4; */
   if ((err = mp_mul_2d(&a2, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 + c; */
   if ((err = mp_add(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a3 << 6; */
   if ((err = mp_mul_2d(&a3, 6, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 + c; */
   if ((err = mp_add(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = a4 << 8; */
   if ((err = mp_mul_2d(&a4, 8, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 + c; */
   if ((err = mp_add(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;

   /** c = b1 << 2; */
   if ((err = mp_mul_2d(&b1, 2, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S6 = b0 + c; */
   if ((err = mp_add(&b0, c, &S6)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b2 << 4; */
   if ((err = mp_mul_2d(&b2, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S6 = S6 + c; */
   if ((err = mp_add(&S6, c, &S6)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b3 << 6; */
   if ((err = mp_mul_2d(&b3, 6, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S6 = S6 + c; */
   if ((err = mp_add(&S6, c, &S6)) != MP_OKAY)                                                           goto LTM_ERR;
   /** c = b4 << 8; */
   if ((err = mp_mul_2d(&b4, 8, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S6 = S6 + c; */
   if ((err = mp_add(&S6, c, &S6)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S4 = S4 * S6; */
   if ((err = mp_mul(&S4, &S6, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S6 = (a0-a1+a2-a3+a4)*(b0-b1+b2-b3+b4); */
   /** S6 = a0 - a1; */
   if ((err = mp_sub(&a0, &a1, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 + a2; */
   if ((err = mp_add(&S6, &a2, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - a3; */
   if ((err = mp_sub(&S6, &a3, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 + a4; */
   if ((err = mp_add(&S6, &a4, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = b0 - b1; */
   if ((err = mp_sub(&b0, &b1, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 + b2; */
   if ((err = mp_add(&S7, &b2, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - b3; */
   if ((err = mp_sub(&S7, &b3, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 + b4; */
   if ((err = mp_add(&S7, &b4, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 * S7; */
   if ((err = mp_mul(&S6, &S7, &S6)) != MP_OKAY)                                                         goto LTM_ERR;

   /** \\S7 = (a0+a1+a2+a3+a4)*(b0+b1+b2+b3+b4); */
   /** S7 = a0 + a1; */
   if ((err = mp_add(&a0, &a1, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 + a2; */
   if ((err = mp_add(&S7, &a2, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 + a3; */
   if ((err = mp_add(&S7, &a3, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 + a4; */
   if ((err = mp_add(&S7, &a4, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = b0 + b1; */
   if ((err = mp_add(&b0, &b1, &S9)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = S9 + b2; */
   if ((err = mp_add(&S9, &b2, &S9)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = S9 + b3; */
   if ((err = mp_add(&S9, &b3, &S9)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = S9 + b4; */
   if ((err = mp_add(&S9, &b4, &S9)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 * S9; */
   if ((err = mp_mul(&S7, &S9, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = a0*b0; */
   if ((err = mp_mul(&a0, &b0, &S9)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - S7; */
   if ((err = mp_sub(&S6, &S7, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 - S5; */
   if ((err = mp_sub(&S2, &S5, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - S9; */
   if ((err = mp_sub(&S4, &S9, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S4 = S4 - (2^16*S1); */
   /** c = S1 << 16; */
   if ((err = mp_mul_2d(&S1, 16, c)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - c; */
   if ((err = mp_sub(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S8 = S8 - S3; */
   if ((err = mp_sub(&S8, &S3, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 / 2; */
   if ((err = mp_div_2(&S6, &S6)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = S5 * 2; */
   if ((err = mp_mul_2(&S5, &S5)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = S5 + S2; */
   if ((err = mp_add(&S5, &S2, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = -S2; */
   if ((err = mp_neg(&S2, &S2)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S8 = -S8; */
   if ((err = mp_neg(&S8, &S8)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S7 = S7 + S6; */
   if ((err = mp_add(&S7, &S6, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = -S6; */
   if ((err = mp_neg(&S6, &S6)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S3 = S3 - S7; */
   if ((err = mp_sub(&S3, &S7, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S5 = S5 - (S7 * 512); */
   /** c = S7 << 9; */
   if ((err = mp_mul_2d(&S7, 9, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S5 = S5 - c; */
   if ((err = mp_sub(&S5, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;

   /** S3 = S3 * 2; */
   if ((err = mp_mul_2(&S3, &S3)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S3 = S3 - S8; */
   if ((err = mp_sub(&S3, &S8, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - S1; */
   if ((err = mp_sub(&S7, &S1, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - S9; */
   if ((err = mp_sub(&S7, &S9, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = S8 + S2; */
   if ((err = mp_add(&S8, &S2, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 + S3; */
   if ((err = mp_add(&S5, &S3, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S8 = S8 - (S6 * 80); */
   /** c = S6 * 80; */
   if ((err = mp_mul_d(&S6, (mp_digit)80, c)) != MP_OKAY)                                                goto LTM_ERR;
   /** S8 = S8 - c; */
   if ((err = mp_sub(&S8, c, &S8)) != MP_OKAY)                                                           goto LTM_ERR;
   /** \\S3 = S3 - (S9 * 510); */
   /** c = S9 * 510; */
   if ((err = mp_mul_d(&S9, (mp_digit)510, c)) != MP_OKAY)                                               goto LTM_ERR;
   /** S3 = S3 - c; */
   if ((err = mp_sub(&S3, c, &S3)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S4 = S4 - S2; */
   if ((err = mp_sub(&S4, &S2, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 * 3; */
   if ((err = mp_mul_d(&S3, (mp_digit)3, &S3)) != MP_OKAY)                                               goto LTM_ERR;
   /** S3 = S3 + S5; */
   if ((err = mp_add(&S3, &S5, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = S8 / 180; */
   if ((err = mp_div_d(&S8, (mp_digit)180, &S8, NULL)) != MP_OKAY)                                       goto LTM_ERR;
   /** \\S5 = S5 + (S7 * 378); */
   /** c = S7 * 378; */
   if ((err = mp_mul_d(&S7, (mp_digit)378, c)) != MP_OKAY)                                               goto LTM_ERR;
   /** S5  = S5 + c; */
   if ((err = mp_add(&S5, c, &S5)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S2 =  S2 >> 2; */
   if ((err = mp_div_2d(&S2, 2, &S2, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S6 = S6 - S2; */
   if ((err = mp_sub(&S6, &S2, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 / (-72); */
   sign = S5.sign;
   if ((err = mp_div_d(&S5, 72, &S5, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   S5.sign = (sign == MP_NEG)?MP_ZPOS: MP_NEG;
   /** S3 = S3 / (-360); */
   sign = S3.sign;
   if ((err = mp_div_d(&S3, (mp_digit)360, &S3, NULL)) != MP_OKAY)                                       goto LTM_ERR;
   S3.sign = (sign == MP_NEG)?MP_ZPOS: MP_NEG;
   /** S2 = S2 - S8; */
   if ((err = mp_sub(&S2, &S8, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - S3; */
   if ((err = mp_sub(&S7, &S3, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S4 = S4 - (S5 * 256); */
   /** c = S5 << 8; */
   if ((err = mp_mul_2d(&S5, 8, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 - c; */
   if ((err = mp_sub(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S3 = S3 - S5; */
   if ((err = mp_sub(&S3, &S5, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S4 = S4 - (S3 * 4096); */
   /** c = S3 << 12; */
   if ((err = mp_mul_2d(&S3, 12, c)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - c; */
   if ((err = mp_sub(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** \\S4 = S4 - (S7 * 16); */
   /** c = S7 << 4; */
   if ((err = mp_mul_2d(&S7, 4, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 - c; */
   if ((err = mp_sub(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** \\S4 = S4 + (S6 * 256); */
   /** c = S6 << 8; */
   if ((err = mp_mul_2d(&S6, 8, c)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 + c; */
   if ((err = mp_add(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S6 = S6 + S2; */
   if ((err = mp_add(&S6, &S2, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 * 180; */
   if ((err = mp_mul_d(&S2, (mp_digit)180, &S2)) != MP_OKAY)                                             goto LTM_ERR;
   /** S2 = S2 + S4; */
   if ((err = mp_add(&S2, &S4, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 / 11340; */
   if ((err = mp_div_d(&S2, (mp_digit)11340, &S2, NULL)) != MP_OKAY)                                     goto LTM_ERR;

   /** \\S4 = S4 + (S6 * 720); */
   /** c = S6 * 720; */
   if ((err = mp_mul_d(&S6, (mp_digit)720, c)) != MP_OKAY)                                               goto LTM_ERR;
   /** S4 = S4 + c; */
   if ((err = mp_add(&S4, c, &S4)) != MP_OKAY)                                                           goto LTM_ERR;
   /** S4 = S4 / (-2160); */
   sign = S4.sign;
   if ((err = mp_div_d(&S4, (mp_digit)2160, &S4, NULL)) != MP_OKAY)                                      goto LTM_ERR;
   S4.sign = (sign == MP_NEG)?MP_ZPOS: MP_NEG;
   /** S6 = S6 - S4; */
   if ((err = mp_sub(&S6, &S4, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = S8 - S2; */
   if ((err = mp_sub(&S8, &S2, &S8)) != MP_OKAY)                                                         goto LTM_ERR;

   /** P = S1*x^8 + S2*x^7 + S3*x^6 + S4*x^5 + S5*x^4 + S6*x^3 + S7*x^2 + S8*x + S9; */
   if ((err = mp_lshd(&S1, 8 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_lshd(&S2, 7 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S3, 6 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S4, 5 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S5, 4 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S5, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S6, 3 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S7, 2 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S7, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S8, 1 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S1, &S8, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S1, &S9, c)) != MP_OKAY)                                                           goto LTM_ERR;
   /** P - A*B \\ == zero */

LTM_ERR:
   mp_clear(&b4);
LTM_ERRb4:
   mp_clear(&b3);
LTM_ERRb3:
   mp_clear(&b2);
LTM_ERRb2:
   mp_clear(&b1);
LTM_ERRb1:
   mp_clear(&b0);
LTM_ERRb0:
   mp_clear(&a4);
LTM_ERRa4:
   mp_clear(&a3);
LTM_ERRa3:
   mp_clear(&a2);
LTM_ERRa2:
   mp_clear(&a1);
LTM_ERRa1:
   mp_clear(&a0);
LTM_ERRa0:
   mp_clear_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, &S8, &S9, NULL);
   return err;
}




#endif
