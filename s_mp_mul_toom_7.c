#include "tommath_private.h"
#ifdef S_MP_MUL_TOOM_7_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */





mp_err s_mp_mul_toom_7(const mp_int *a, const mp_int *b,  mp_int *c)
{
   mp_int S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, t;
   mp_int a0, a1, a2, a3, a4, a5, a6;
   mp_int b0, b1, b2, b3, b4, b5, b6;
   mp_err err = MP_OKAY;
   int B;

   B = a->used / 7;

   if ((err = mp_init_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, &S8, &S9, &S10, &S11, &S12, &S13, &t, NULL)) != MP_OKAY) {
      return err;
   }

   if ((err = mp_init_size(&a0, B)) != MP_OKAY)                                                          goto LTM_ERRa0;
   if ((err = mp_init_size(&a1, B)) != MP_OKAY)                                                          goto LTM_ERRa1;
   if ((err = mp_init_size(&a2, B)) != MP_OKAY)                                                          goto LTM_ERRa2;
   if ((err = mp_init_size(&a3, B)) != MP_OKAY)                                                          goto LTM_ERRa3;
   if ((err = mp_init_size(&a4, B)) != MP_OKAY)                                                          goto LTM_ERRa4;
   if ((err = mp_init_size(&a5, B)) != MP_OKAY)                                                          goto LTM_ERRa5;
   if ((err = mp_init_size(&a6, a->used - 6 * B)) != MP_OKAY)                                            goto LTM_ERRa6;


   if ((err = mp_init_size(&b0, B)) != MP_OKAY)                                                          goto LTM_ERRb0;
   if ((err = mp_init_size(&b1, B)) != MP_OKAY)                                                          goto LTM_ERRb1;
   if ((err = mp_init_size(&b2, B)) != MP_OKAY)                                                          goto LTM_ERRb2;
   if ((err = mp_init_size(&b3, B)) != MP_OKAY)                                                          goto LTM_ERRb3;
   if ((err = mp_init_size(&b4, B)) != MP_OKAY)                                                          goto LTM_ERRb4;
   if ((err = mp_init_size(&b5, B)) != MP_OKAY)                                                          goto LTM_ERRb5;
   if ((err = mp_init_size(&b6, b->used - 6 * B)) != MP_OKAY)                                            goto LTM_ERRb6;

   /** A = a0 + a1 * x + a2 * x^2 + a3 * x^3 + a4 * x^4 + a5 * x^5 + a6 * x^6; */
   a0.used = a1.used = a2.used = a3.used = a4.used = a5.used = B;
   a6.used = a->used - 6 * B;
   s_mp_copy_digs(a0.dp, a->dp, a0.used);
   s_mp_copy_digs(a1.dp, a->dp + B, a1.used);
   s_mp_copy_digs(a2.dp, a->dp + 2 * B, a2.used);
   s_mp_copy_digs(a3.dp, a->dp + 3 * B, a3.used);
   s_mp_copy_digs(a4.dp, a->dp + 4 * B, a4.used);
   s_mp_copy_digs(a5.dp, a->dp + 5 * B, a5.used);
   s_mp_copy_digs(a6.dp, a->dp + 6 * B, a6.used);
   mp_clamp(&a0);
   mp_clamp(&a1);
   mp_clamp(&a2);
   mp_clamp(&a3);
   mp_clamp(&a4);
   mp_clamp(&a5);
   mp_clamp(&a6);

   /** B = b0 + b1 * x + b2 * x^2 + b3 * x^3 + b4 * x^4 + b5 * x^5 + b6 * x^6; */
   b0.used = b1.used = b2.used = b3.used = b4.used = b5.used = B;
   b6.used = b->used - 6 * B;
   s_mp_copy_digs(b0.dp, b->dp, b0.used);
   s_mp_copy_digs(b1.dp, b->dp + B, b1.used);
   s_mp_copy_digs(b2.dp, b->dp + 2 * B, b2.used);
   s_mp_copy_digs(b3.dp, b->dp + 3 * B, b3.used);
   s_mp_copy_digs(b4.dp, b->dp + 4 * B, b4.used);
   s_mp_copy_digs(b5.dp, b->dp + 5 * B, b5.used);
   s_mp_copy_digs(b6.dp, b->dp + 6 * B, b6.used);
   mp_clamp(&b0);
   mp_clamp(&b1);
   mp_clamp(&b2);
   mp_clamp(&b3);
   mp_clamp(&b4);
   mp_clamp(&b5);
   mp_clamp(&b6);


   /* Evaluation */

   /** S1 = a0 + a2; */
   if ((err = mp_add(&a0, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + a4; */
   if ((err = mp_add(&S1, &a4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + a6; */
   if ((err = mp_add(&S1, &a6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = a1 + a3; */
   if ((err = mp_add(&a1, &a3, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S13 = S13 + a5; */
   if ((err = mp_add(&S13, &a5, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S4 = S1 - S13; */
   if ((err = mp_sub(&S1, &S13, &S4)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S5 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S5)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = b0 + b2; */
   if ((err = mp_add(&b0, &b2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + b4; */
   if ((err = mp_add(&S1, &b4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + b6; */
   if ((err = mp_add(&S1, &b6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = b1 + b3; */
   if ((err = mp_add(&b1, &b3, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S13 = S13 + b5; */
   if ((err = mp_add(&S13, &b5, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S6 = S1 - S13; */
   if ((err = mp_sub(&S1, &S13, &S6)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S7)) != MP_OKAY)                                                        goto LTM_ERR;

   /** S4 = S4 * S6; */
   if ((err = mp_mul(&S4, &S6, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 * S7; */
   if ((err = mp_mul(&S5, &S7, &S5)) != MP_OKAY)                                                         goto LTM_ERR;

   /** S1 = a0<<2; */
   if ((err = mp_mul_2d(&a0, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a2; */
   if ((err = mp_add(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<2; */
   if ((err = mp_mul_2d(&S1, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a4; */
   if ((err = mp_add(&S1, &a4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<2; */
   if ((err = mp_mul_2d(&S1, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a6; */
   if ((err = mp_add(&S1, &a6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = a1<<2; */
   if ((err = mp_mul_2d(&a1, 2, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 + a3; */
   if ((err = mp_add(&S13, &a3, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<2; */
   if ((err = mp_mul_2d(&S13, 2, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 + a5; */
   if ((err = mp_add(&S13, &a5, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<1; */
   if ((err = mp_mul_2d(&S13, 1, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S3 = S1 - S13; */
   if ((err = mp_sub(&S1, &S13, &S3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S6 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S6)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = b0<<2; */
   if ((err = mp_mul_2d(&b0, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b2; */
   if ((err = mp_add(&S1, &b2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<2; */
   if ((err = mp_mul_2d(&S1, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b4; */
   if ((err = mp_add(&S1, &b4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<2; */
   if ((err = mp_mul_2d(&S1, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b6; */
   if ((err = mp_add(&S1, &b6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = b1<<2; */
   if ((err = mp_mul_2d(&b1, 2, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 + b3; */
   if ((err = mp_add(&S13, &b3, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<2; */
   if ((err = mp_mul_2d(&S13, 2, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 + b5; */
   if ((err = mp_add(&S13, &b5, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<1; */
   if ((err = mp_mul_2d(&S13, 1, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S7 = S1 - S13; */
   if ((err = mp_sub(&S1, &S13, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S8 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S8)) != MP_OKAY)                                                        goto LTM_ERR;

   /** S3 = S3 * S7; */
   if ((err = mp_mul(&S3, &S7, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 * S8; */
   if ((err = mp_mul(&S6, &S8, &S6)) != MP_OKAY)                                                         goto LTM_ERR;

   /** S1 = a0<<4; */
   if ((err = mp_mul_2d(&a0, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a2; */
   if ((err = mp_add(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a4; */
   if ((err = mp_add(&S1, &a4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a6; */
   if ((err = mp_add(&S1, &a6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = a1<<4; */
   if ((err = mp_mul_2d(&a1, 4, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 + a3; */
   if ((err = mp_add(&S13, &a3, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<4; */
   if ((err = mp_mul_2d(&S13, 4, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 + a5; */
   if ((err = mp_add(&S13, &a5, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<2; */
   if ((err = mp_mul_2d(&S13, 2, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S7 = S1 - S13; */
   if ((err = mp_sub(&S1, &S13, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S8 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S8)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = b0<<4; */
   if ((err = mp_mul_2d(&b0, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b2; */
   if ((err = mp_add(&S1, &b2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b4; */
   if ((err = mp_add(&S1, &b4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b6; */
   if ((err = mp_add(&S1, &b6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = b1<<4; */
   if ((err = mp_mul_2d(&b1, 4, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 + b3; */
   if ((err = mp_add(&S13, &b3, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<4; */
   if ((err = mp_mul_2d(&S13, 4, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 + b5; */
   if ((err = mp_add(&S13, &b5, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<2; */
   if ((err = mp_mul_2d(&S13, 2, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S9 = S1 - S13; */
   if ((err = mp_sub(&S1, &S13, &S9)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S10 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S10)) != MP_OKAY)                                                       goto LTM_ERR;

   /** S2 = S7 * S9; */
   if ((err = mp_mul(&S7, &S9, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S8 * S10; */
   if ((err = mp_mul(&S8, &S10, &S7)) != MP_OKAY)                                                        goto LTM_ERR;

   /** S1 = a6<<2; */
   if ((err = mp_mul_2d(&a6, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a4; */
   if ((err = mp_add(&S1, &a4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<2; */
   if ((err = mp_mul_2d(&S1, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a2; */
   if ((err = mp_add(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<2; */
   if ((err = mp_mul_2d(&S1, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a0; */
   if ((err = mp_add(&S1, &a0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = a5<<2; */
   if ((err = mp_mul_2d(&a5, 2, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 + a3; */
   if ((err = mp_add(&S13, &a3, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<2; */
   if ((err = mp_mul_2d(&S13, 2, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 + a1; */
   if ((err = mp_add(&S13, &a1, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<1; */
   if ((err = mp_mul_2d(&S13, 1, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S8 = S1 - S13; */
   if ((err = mp_sub(&S1, &S13, &S8)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S9 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S9)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = b6<<2; */
   if ((err = mp_mul_2d(&b6, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b4; */
   if ((err = mp_add(&S1, &b4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<2; */
   if ((err = mp_mul_2d(&S1, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b2; */
   if ((err = mp_add(&S1, &b2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<2; */
   if ((err = mp_mul_2d(&S1, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b0; */
   if ((err = mp_add(&S1, &b0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = b5<<2; */
   if ((err = mp_mul_2d(&b5, 2, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 + b3; */
   if ((err = mp_add(&S13, &b3, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<2; */
   if ((err = mp_mul_2d(&S13, 2, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 + b1; */
   if ((err = mp_add(&S13, &b1, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<1; */
   if ((err = mp_mul_2d(&S13, 1, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S10 = S1 - S13; */
   if ((err = mp_sub(&S1, &S13, &S10)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S11 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S11)) != MP_OKAY)                                                       goto LTM_ERR;

   /** S8 = S8 * S10; */
   if ((err = mp_mul(&S8, &S10, &S8)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S9 = S9 * S11; */
   if ((err = mp_mul(&S9, &S11, &S9)) != MP_OKAY)                                                        goto LTM_ERR;

   /** S1 = a6<<4; */
   if ((err = mp_mul_2d(&a6, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a4; */
   if ((err = mp_add(&S1, &a4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a2; */
   if ((err = mp_add(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a0; */
   if ((err = mp_add(&S1, &a0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = a5<<4; */
   if ((err = mp_mul_2d(&a5, 4, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 + a3; */
   if ((err = mp_add(&S13, &a3, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<4; */
   if ((err = mp_mul_2d(&S13, 4, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 + a1; */
   if ((err = mp_add(&S13, &a1, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<2; */
   if ((err = mp_mul_2d(&S13, 2, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S10 = S1 - S13; */
   if ((err = mp_sub(&S1, &S13, &S10)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S11 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S11)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S1 = b6<<4; */
   if ((err = mp_mul_2d(&b6, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b4; */
   if ((err = mp_add(&S1, &b4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b2; */
   if ((err = mp_add(&S1, &b2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b0; */
   if ((err = mp_add(&S1, &b0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = b5<<4; */
   if ((err = mp_mul_2d(&b5, 4, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 + b3; */
   if ((err = mp_add(&S13, &b3, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<4; */
   if ((err = mp_mul_2d(&S13, 4, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 + b1; */
   if ((err = mp_add(&S13, &b1, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<2; */
   if ((err = mp_mul_2d(&S13, 2, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S12 = S1 - S13; */
   if ((err = mp_sub(&S1, &S13, &S12)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S13)) != MP_OKAY)                                                       goto LTM_ERR;

   /** S10 = S10 * S12; */
   if ((err = mp_mul(&S10, &S12, &S10)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S11 = S11 * S13; */
   if ((err = mp_mul(&S11, &S13, &S11)) != MP_OKAY)                                                      goto LTM_ERR;

   /** S1 = a6<<6; */
   if ((err = mp_mul_2d(&a6, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a4; */
   if ((err = mp_add(&S1, &a4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<6; */
   if ((err = mp_mul_2d(&S1, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a2; */
   if ((err = mp_add(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<6; */
   if ((err = mp_mul_2d(&S1, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a0; */
   if ((err = mp_add(&S1, &a0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = a5<<6; */
   if ((err = mp_mul_2d(&a5, 6, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 + a3; */
   if ((err = mp_add(&S13, &a3, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<6; */
   if ((err = mp_mul_2d(&S13, 6, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 + a1; */
   if ((err = mp_add(&S13, &a1, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<3; */
   if ((err = mp_mul_2d(&S13, 3, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S12 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S12)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S1 = b6<<6; */
   if ((err = mp_mul_2d(&b6, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b4; */
   if ((err = mp_add(&S1, &b4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<6; */
   if ((err = mp_mul_2d(&S1, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b2; */
   if ((err = mp_add(&S1, &b2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<6; */
   if ((err = mp_mul_2d(&S1, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + b0; */
   if ((err = mp_add(&S1, &b0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = b5<<6; */
   if ((err = mp_mul_2d(&b5, 6, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 + b3; */
   if ((err = mp_add(&S13, &b3, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<6; */
   if ((err = mp_mul_2d(&S13, 6, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 + b1; */
   if ((err = mp_add(&S13, &b1, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13<<3; */
   if ((err = mp_mul_2d(&S13, 3, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S1 + S13; */
   if ((err = mp_add(&S1, &S13, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S12 = S12 * S13; */
   if ((err = mp_mul(&S12, &S13, &S12)) != MP_OKAY)                                                      goto LTM_ERR;

   /** S1 = a0 * b0; */
   if ((err = mp_mul(&a0, &b0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S13 = a6 * b6; */
   if ((err = mp_mul(&a6, &b6, &S13)) != MP_OKAY)                                                        goto LTM_ERR;

   /* Interpolation */

   /** S2 = -S2; */
   S2.sign = (S2.sign == MP_NEG)?MP_ZPOS: MP_NEG;
   /** S2 = S2 + S7; */
   if ((err = mp_add(&S2, &S7, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S1<<24; */
   if ((err = mp_mul_2d(&S1, 24, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = S7 - t; */
   if ((err = mp_sub(&S7, &t, &S7)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S2 = S2 / 2; */
   if ((err = mp_div_2(&S2, &S2)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S7 = S7 - S2; */
   if ((err = mp_sub(&S7, &S2, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - S13; */
   if ((err = mp_sub(&S7, &S13, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S3 - S6; */
   if ((err = mp_sub(&S3, &S6, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = -S3; */
   S3.sign = (S3.sign == MP_NEG)?MP_ZPOS: MP_NEG;
   /** S3 = S3 / 2; */
   if ((err = mp_div_2(&S3, &S3)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S6 = S6 - S3; */
   if ((err = mp_sub(&S6, &S3, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - S13; */
   if ((err = mp_sub(&S6, &S13, &S6)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S1<<12; */
   if ((err = mp_mul_2d(&S1, 12, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S6 = S6 - t; */
   if ((err = mp_sub(&S6, &t, &S6)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S11 = S11 - S10; */
   if ((err = mp_sub(&S11, &S10, &S11)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S11 = S11 / 2; */
   if ((err = mp_div_2(&S11, &S11)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S10 = S10 + S11; */
   if ((err = mp_add(&S10, &S11, &S10)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S10 = S10 - S1; */
   if ((err = mp_sub(&S10, &S1, &S10)) != MP_OKAY)                                                       goto LTM_ERR;
   /** t = S13<<24; */
   if ((err = mp_mul_2d(&S13, 24, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S9 = S9 - S8; */
   if ((err = mp_sub(&S9, &S8, &S9)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = S9 / 2; */
   if ((err = mp_div_2(&S9, &S9)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S8 = S8 + S9; */
   if ((err = mp_add(&S8, &S9, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = S8 - S1; */
   if ((err = mp_sub(&S8, &S1, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S13<<12; */
   if ((err = mp_mul_2d(&S13, 12, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S8 = S8 - t; */
   if ((err = mp_sub(&S8, &t, &S8)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S13<<36; */
   if ((err = mp_mul_2d(&S13, 36, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S12 = S12 - S1; */
   if ((err = mp_sub(&S12, &S1, &S12)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S4 = S4 + S5; */
   if ((err = mp_add(&S4, &S5, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 / 2; */
   if ((err = mp_div_2(&S4, &S4)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S5 = S5 - S4; */
   if ((err = mp_sub(&S5, &S4, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - S1; */
   if ((err = mp_sub(&S4, &S1, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - S13; */
   if ((err = mp_sub(&S4, &S13, &S4)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S2 = S2 + S11; */
   if ((err = mp_add(&S2, &S11, &S2)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = S7 + S10; */
   if ((err = mp_add(&S7, &S10, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S3 + S9; */
   if ((err = mp_add(&S3, &S9, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 + S8; */
   if ((err = mp_add(&S6, &S8, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S4<<13; */
   if ((err = mp_mul_2d(&S4, 13, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = S7 - t; */
   if ((err = mp_sub(&S7, &t, &S7)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S4<<7; */
   if ((err = mp_mul_2d(&S4, 7, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - t; */
   if ((err = mp_sub(&S6, &t, &S6)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S5 * 17408; */
   if ((err = mp_mul_d(&S5, 17408, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S5 * 160; */
   if ((err = mp_mul_d(&S5, 160, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S3 - t; */
   if ((err = mp_sub(&S3, &t, &S3)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S6 * 400; */
   if ((err = mp_mul_d(&S6, 400, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = S7 - t; */
   if ((err = mp_sub(&S7, &t, &S7)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S3 * 680; */
   if ((err = mp_mul_d(&S3, 680, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S2 = S2 / 2891700; */
   if ((err = mp_div_d(&S2, 2891700, &S2, NULL)) != MP_OKAY)                                             goto LTM_ERR;
   /** S7 = S7 / 680400; */
   if ((err = mp_div_d(&S7, 680400, &S7, NULL)) != MP_OKAY)                                              goto LTM_ERR;
   /** t = S2 * 1890; */
   if ((err = mp_mul_d(&S2, 1890, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S3 = S3 - t; */
   if ((err = mp_sub(&S3, &t, &S3)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = S3 / 360; */
   if ((err = mp_div_d(&S3, 360, &S3, NULL)) != MP_OKAY)                                                 goto LTM_ERR;
   /** t = S7 * 900; */
   if ((err = mp_mul_d(&S7, 900, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S6 = S6 - t; */
   if ((err = mp_sub(&S6, &t, &S6)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S6 = S6 / 144; */
   if ((err = mp_div_d(&S6, 144, &S6, NULL)) != MP_OKAY)                                                 goto LTM_ERR;
   /** S5 = S5 - S2; */
   if ((err = mp_sub(&S5, &S2, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 - S3; */
   if ((err = mp_sub(&S5, &S3, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - S7; */
   if ((err = mp_sub(&S4, &S7, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - S6; */
   if ((err = mp_sub(&S4, &S6, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S4<<6; */
   if ((err = mp_mul_2d(&S4, 6, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = S8 - t; */
   if ((err = mp_sub(&S8, &t, &S8)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S4<<12; */
   if ((err = mp_mul_2d(&S4, 12, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S4<<18; */
   if ((err = mp_mul_2d(&S4, 18, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S5 * 32; */
   if ((err = mp_mul_d(&S5, 32, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S5<<10; */
   if ((err = mp_mul_2d(&S5, 10, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S5<<15; */
   if ((err = mp_mul_2d(&S5, 15, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S6<<4; */
   if ((err = mp_mul_2d(&S6, 4, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = S8 - t; */
   if ((err = mp_sub(&S8, &t, &S8)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S6<<8; */
   if ((err = mp_mul_2d(&S6, 8, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S6<<12; */
   if ((err = mp_mul_2d(&S6, 12, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S3<<3; */
   if ((err = mp_mul_2d(&S3, 3, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S3<<6; */
   if ((err = mp_mul_2d(&S3, 6, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S3<<9; */
   if ((err = mp_mul_2d(&S3, 9, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S7<<6; */
   if ((err = mp_mul_2d(&S7, 6, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S7<<2; */
   if ((err = mp_mul_2d(&S7, 2, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = S8 - t; */
   if ((err = mp_sub(&S8, &t, &S8)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S7<<4; */
   if ((err = mp_mul_2d(&S7, 4, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S2 * 2; */
   if ((err = mp_mul_2(&S2, &t)) != MP_OKAY)                                                             goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S2 * 4; */
   if ((err = mp_mul_d(&S2, 4, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S2<<3; */
   if ((err = mp_mul_2d(&S2, 3, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S8 * 272; */
   if ((err = mp_mul_d(&S8, 272, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S10 = S10 / 771120; */
   if ((err = mp_div_d(&S10, 771120, &S10, NULL)) != MP_OKAY)                                            goto LTM_ERR;
   /** t = S10 * 1020; */
   if ((err = mp_mul_d(&S10, 1020, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S8 = S8 - t; */
   if ((err = mp_sub(&S8, &t, &S8)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S8 = S8 / 240; */
   if ((err = mp_div_d(&S8, 240, &S8, NULL)) != MP_OKAY)                                                 goto LTM_ERR;
   /** t = S8 * 16773120; */
   if ((err = mp_mul_d(&S8, 16773120, &t)) != MP_OKAY)                                                   goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S10<<6; */
   if ((err = mp_mul_2d(&S10, 6, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t * 16777215; */
   if ((err = mp_mul_d(&t, 16777215, &t)) != MP_OKAY)                                                    goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S9 * 21504; */
   if ((err = mp_mul_d(&S9, 21504, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S9 * 160; */
   if ((err = mp_mul_d(&S9, 160, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S11 * 2210; */
   if ((err = mp_mul_d(&S11, 2210, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S12 = -S12; */
   S12.sign = (S12.sign == MP_NEG)?MP_ZPOS: MP_NEG;
   /** S12 = S12>>7; */
   if ((err = mp_div_2d(&S12, 7, &S12, NULL)) != MP_OKAY)                                                goto LTM_ERR;
   /** S12 = S12 / 2168775; */
   if ((err = mp_div_d(&S12, 2168775, &S12, NULL)) != MP_OKAY)                                           goto LTM_ERR;
   /** t = S12 * 181440; */
   if ((err = mp_mul_d(&S12, 181440, &t)) != MP_OKAY)                                                    goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S11 = S11 / 3866940; */
   if ((err = mp_div_d(&S11, 3866940, &S11, NULL)) != MP_OKAY)                                           goto LTM_ERR;
   /** t = S12 * 504; */
   if ((err = mp_mul_d(&S12, 504, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S11 * 2046; */
   if ((err = mp_mul_d(&S11, 2046, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S9 = S9 / 96; */
   if ((err = mp_div_d(&S9, 96, &S9, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S2 = S2 - S11; */
   if ((err = mp_sub(&S2, &S11, &S2)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S3 - S12; */
   if ((err = mp_sub(&S3, &S12, &S3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S5 = S5 - S9; */
   if ((err = mp_sub(&S5, &S9, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - S8; */
   if ((err = mp_sub(&S6, &S8, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - S10; */
   if ((err = mp_sub(&S7, &S10, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S3;S3 = S7;S7 = t;     */
   mp_exch(&S3, &S7);
   /** t = S4;S4 = S7;S7 = t;     */
   mp_exch(&S4, &S7);
   /** t = S5;S5 = S6;S6 = t;     */
   mp_exch(&S5, &S6);
   /** t = S8;S8 = S9;S9 = t;     */
   mp_exch(&S8, &S9);
   /** t = S10;S10 = S12;S12 = t; */
   mp_exch(&S10, &S12);
   /** t = S11;S11 = S12;S12 = t; */
   mp_exch(&S11, &S12);




   /** P =  S1 * x^0+  S2 * x^1 + S3 * x^2 + S4  * x^3 + S5  * x^4  + S6  * x^5; */
   /** P += S7 * x^6 + S8 * x^7 + S9 * x^8 + S10 * x^9 + S11 * x^10 + S12 * x^11 + S13 * x^12; */


   if ((err = mp_lshd(&S2, 1 *  B)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_add(&S1, &S2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S3, 2 *  B)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_add(&S1, &S3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S4, 3 *  B)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_add(&S1, &S4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S5, 4 *  B)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_add(&S1, &S5, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S6, 5 *  B)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_add(&S1, &S6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S7, 6 *  B)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_add(&S1, &S7, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S8, 7 *  B)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_add(&S1, &S8, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S9, 8 *  B)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_add(&S1, &S9, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S10, 9 *  B)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S1, &S10, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_lshd(&S11, 10 *  B)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_add(&S1, &S11, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_lshd(&S12, 11 *  B)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_add(&S1, &S12, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_lshd(&S13, 12 *  B)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_add(&S1, &S13, &S1)) != MP_OKAY)                                                        goto LTM_ERR;

   mp_exch(&S1, c);

   /** P - A*B */

LTM_ERR:
   mp_clear(&b6);
LTM_ERRb6:
   mp_clear(&b5);
LTM_ERRb5:
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
   mp_clear(&a6);
LTM_ERRa6:
   mp_clear(&a5);
LTM_ERRa5:
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
   mp_clear_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, &S8, &S9, &S10, &S11, &S12, &S13, &t, NULL);
   return err;
}

#endif
