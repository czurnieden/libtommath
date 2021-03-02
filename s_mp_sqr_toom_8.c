#include "tommath_private.h"
#ifdef S_MP_SQR_TOOM_8_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */


mp_err s_mp_sqr_toom_8(const mp_int *a, mp_int *c)
{
   mp_int S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, t;
   mp_int a0, a1, a2, a3, a4, a5, a6, a7;
   mp_err err = MP_OKAY;
   int B;

   B = a->used / 8;

   if ((err = mp_init_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, &S8, &S9, &S10, &S11, &S12, &S13, &S14, &S15, &t,
                            NULL)) != MP_OKAY) {
      return err;
   }

   if ((err = mp_init_size(&a0, B)) != MP_OKAY)                                                          goto LTM_ERRa0;
   if ((err = mp_init_size(&a1, B)) != MP_OKAY)                                                          goto LTM_ERRa1;
   if ((err = mp_init_size(&a2, B)) != MP_OKAY)                                                          goto LTM_ERRa2;
   if ((err = mp_init_size(&a3, B)) != MP_OKAY)                                                          goto LTM_ERRa3;
   if ((err = mp_init_size(&a4, B)) != MP_OKAY)                                                          goto LTM_ERRa4;
   if ((err = mp_init_size(&a5, B)) != MP_OKAY)                                                          goto LTM_ERRa5;
   if ((err = mp_init_size(&a6, B)) != MP_OKAY)                                                          goto LTM_ERRa6;
   if ((err = mp_init_size(&a7, a->used - 7 * B)) != MP_OKAY)                                            goto LTM_ERRa7;

   /* A = a0 * x^0 + a1 * x^1 + a2 * x^2 + a3 * x^3  + a4 * x^4 + a5 * x^5 + a6 * x^6 + a7 * x^7; */
   a0.used = a1.used = a2.used = a3.used = a4.used = a5.used = a6.used = B;
   a7.used = a->used - 7 * B;
   s_mp_copy_digs(a0.dp, a->dp, a0.used);
   s_mp_copy_digs(a1.dp, a->dp + B, a1.used);
   s_mp_copy_digs(a2.dp, a->dp + 2 * B, a2.used);
   s_mp_copy_digs(a3.dp, a->dp + 3 * B, a3.used);
   s_mp_copy_digs(a4.dp, a->dp + 4 * B, a4.used);
   s_mp_copy_digs(a5.dp, a->dp + 5 * B, a5.used);
   s_mp_copy_digs(a6.dp, a->dp + 6 * B, a6.used);
   s_mp_copy_digs(a7.dp, a->dp + 7 * B, a7.used);
   mp_clamp(&a0);
   mp_clamp(&a1);
   mp_clamp(&a2);
   mp_clamp(&a3);
   mp_clamp(&a4);
   mp_clamp(&a5);
   mp_clamp(&a6);
   mp_clamp(&a7);


   /** S1 = a0 + a2; */
   if ((err = mp_add(&a0, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + a4; */
   if ((err = mp_add(&S1, &a4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + a6; */
   if ((err = mp_add(&S1, &a6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S15 = a1 + a3; */
   if ((err = mp_add(&a1, &a3, &S15)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S15 = S15 + a5; */
   if ((err = mp_add(&S15, &a5, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15 + a7; */
   if ((err = mp_add(&S15, &a7, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S2 = S1 - S15; */
   if ((err = mp_sub(&S1, &S15, &S2)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S1 + S15; */
   if ((err = mp_add(&S1, &S15, &S3)) != MP_OKAY)                                                        goto LTM_ERR;

   /** S5 = S2^2; */
   if ((err = mp_sqr(&S2, &S5)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S6 = S3^2; */
   if ((err = mp_sqr(&S3, &S6)) != MP_OKAY)                                                              goto LTM_ERR;

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
   /** S1 = S1<<2; */
   if ((err = mp_mul_2d(&S1, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S15 = a1<<4; */
   if ((err = mp_mul_2d(&a1, 4, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15 + a3; */
   if ((err = mp_add(&S15, &a3, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<4; */
   if ((err = mp_mul_2d(&S15, 4, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a5; */
   if ((err = mp_add(&S15, &a5, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<4; */
   if ((err = mp_mul_2d(&S15, 4, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a7; */
   if ((err = mp_add(&S15, &a7, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S11 = S1 - S15; */
   if ((err = mp_sub(&S1, &S15, &S11)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S12 = S1 + S15; */
   if ((err = mp_add(&S1, &S15, &S12)) != MP_OKAY)                                                       goto LTM_ERR;

   /** S3 = S11^2; */
   if ((err = mp_sqr(&S11, &S3)) != MP_OKAY)                                                             goto LTM_ERR;
   /** S8 = S12^2; */
   if ((err = mp_sqr(&S12, &S8)) != MP_OKAY)                                                             goto LTM_ERR;

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
   /** S1 = S1<<1; */
   if ((err = mp_mul_2d(&S1, 1, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S15 = a1<<2; */
   if ((err = mp_mul_2d(&a1, 2, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15 + a3; */
   if ((err = mp_add(&S15, &a3, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<2; */
   if ((err = mp_mul_2d(&S15, 2, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a5; */
   if ((err = mp_add(&S15, &a5, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<2; */
   if ((err = mp_mul_2d(&S15, 2, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a7; */
   if ((err = mp_add(&S15, &a7, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S11 = S1 - S15; */
   if ((err = mp_sub(&S1, &S15, &S11)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S12 = S1 + S15; */
   if ((err = mp_add(&S1, &S15, &S12)) != MP_OKAY)                                                       goto LTM_ERR;

   /** S4 = S11^2; */
   if ((err = mp_sqr(&S11, &S4)) != MP_OKAY)                                                             goto LTM_ERR;
   /** S7 = S12^2; */
   if ((err = mp_sqr(&S12, &S7)) != MP_OKAY)                                                             goto LTM_ERR;

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
   /** S15 = a7<<2; */
   if ((err = mp_mul_2d(&a7, 2, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15 + a5; */
   if ((err = mp_add(&S15, &a5, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<2; */
   if ((err = mp_mul_2d(&S15, 2, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a3; */
   if ((err = mp_add(&S15, &a3, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<2; */
   if ((err = mp_mul_2d(&S15, 2, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a1; */
   if ((err = mp_add(&S15, &a1, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<1; */
   if ((err = mp_mul_2d(&S15, 1, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S11 = S1 + S15; */
   if ((err = mp_add(&S1, &S15, &S11)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S12 = S1 - S15; */
   if ((err = mp_sub(&S1, &S15, &S12)) != MP_OKAY)                                                       goto LTM_ERR;

   /** S9 = S12^2; */
   if ((err = mp_sqr(&S12, &S9)) != MP_OKAY)                                                             goto LTM_ERR;
   /** S10 = S11^2; */
   if ((err = mp_sqr(&S11, &S10)) != MP_OKAY)                                                            goto LTM_ERR;

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
   /** S15 = a7<<4; */
   if ((err = mp_mul_2d(&a7, 4, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15 + a5; */
   if ((err = mp_add(&S15, &a5, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<4; */
   if ((err = mp_mul_2d(&S15, 4, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a3; */
   if ((err = mp_add(&S15, &a3, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<4; */
   if ((err = mp_mul_2d(&S15, 4, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a1; */
   if ((err = mp_add(&S15, &a1, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<2; */
   if ((err = mp_mul_2d(&S15, 2, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S1 + S15; */
   if ((err = mp_add(&S1, &S15, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S14 = S1 - S15; */
   if ((err = mp_sub(&S1, &S15, &S14)) != MP_OKAY)                                                       goto LTM_ERR;

   /** S12 = S13^2; */
   if ((err = mp_sqr(&S13, &S12)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S11 = S14^2; */
   if ((err = mp_sqr(&S14, &S11)) != MP_OKAY)                                                            goto LTM_ERR;

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
   /** S15 = a7<<6; */
   if ((err = mp_mul_2d(&a7, 6, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15 + a5; */
   if ((err = mp_add(&S15, &a5, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<6; */
   if ((err = mp_mul_2d(&S15, 6, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a3; */
   if ((err = mp_add(&S15, &a3, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<6; */
   if ((err = mp_mul_2d(&S15, 6, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a1; */
   if ((err = mp_add(&S15, &a1, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<3; */
   if ((err = mp_mul_2d(&S15, 3, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S14 = S1 + S15; */
   if ((err = mp_add(&S1, &S15, &S14)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S1 - S15; */
   if ((err = mp_sub(&S1, &S15, &S13)) != MP_OKAY)                                                       goto LTM_ERR;


   /** S14 = S14^2; */
   if ((err = mp_sqr(&S14, &S14)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S13 = S13^2; */
   if ((err = mp_sqr(&S13, &S13)) != MP_OKAY)                                                            goto LTM_ERR;

   /** S1 = a0<<6; */
   if ((err = mp_mul_2d(&a0, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a2; */
   if ((err = mp_add(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<6; */
   if ((err = mp_mul_2d(&S1, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a4; */
   if ((err = mp_add(&S1, &a4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<6; */
   if ((err = mp_mul_2d(&S1, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a6; */
   if ((err = mp_add(&S1, &a6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<3; */
   if ((err = mp_mul_2d(&S1, 3, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S15 = a1<<6; */
   if ((err = mp_mul_2d(&a1, 6, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15 + a3; */
   if ((err = mp_add(&S15, &a3, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<6; */
   if ((err = mp_mul_2d(&S15, 6, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a5; */
   if ((err = mp_add(&S15, &a5, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S15<<6; */
   if ((err = mp_mul_2d(&S15, 6, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 + a7; */
   if ((err = mp_add(&S15, &a7, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S2 = S1 - S15; */
   if ((err = mp_sub(&S1, &S15, &S2)) != MP_OKAY)                                                        goto LTM_ERR;

   /** S2 = S2^2; */
   if ((err = mp_sqr(&S2, &S2)) != MP_OKAY)                                                              goto LTM_ERR;

   /** S1 = a0^2; */
   if ((err = mp_sqr(&a0, &S1)) != MP_OKAY)                                                              goto LTM_ERR;

   /** S15 = a7^2; */
   if ((err = mp_sqr(&a7, &S15)) != MP_OKAY)                                                             goto LTM_ERR;


   /* TODO: interpolation is the same for 8-way multiplication and 8-way square  */
   /** S3 = S3 + S8; */
   if ((err = mp_add(&S3, &S8, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3>>1; */
   if ((err = mp_div_2d(&S3, 1, &S3, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S8 = S8 - S3; */
   if ((err = mp_sub(&S8, &S3, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 + S7; */
   if ((err = mp_add(&S4, &S7, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4>>1; */
   if ((err = mp_div_2d(&S4, 1, &S4, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S7 = S7 - S4; */
   if ((err = mp_sub(&S7, &S4, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = S9 + S10; */
   if ((err = mp_add(&S9, &S10, &S9)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S9 = S9>>1; */
   if ((err = mp_div_2d(&S9, 1, &S9, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S10 = S10 - S9; */
   if ((err = mp_sub(&S10, &S9, &S10)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S11 = S11 + S12; */
   if ((err = mp_add(&S11, &S12, &S11)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S11 = S11>>1; */
   if ((err = mp_div_2d(&S11, 1, &S11, NULL)) != MP_OKAY)                                                goto LTM_ERR;
   /** S12 = S12 - S11; */
   if ((err = mp_sub(&S12, &S11, &S12)) != MP_OKAY)                                                      goto LTM_ERR;
   /** t = S1<<42; */
   if ((err = mp_mul_2d(&S1, 42, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S15<<42; */
   if ((err = mp_mul_2d(&S15, 42, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S2 = S2 + S13; */
   if ((err = mp_add(&S2, &S13, &S2)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S13 = S13 + S14; */
   if ((err = mp_add(&S13, &S14, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13>>1; */
   if ((err = mp_div_2d(&S13, 1, &S13, NULL)) != MP_OKAY)                                                goto LTM_ERR;
   /** S14 = S14 - S13; */
   if ((err = mp_sub(&S14, &S13, &S14)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S5 = S5 + S6; */
   if ((err = mp_add(&S5, &S6, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5>>1; */
   if ((err = mp_div_2d(&S5, 1, &S5, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S6 = S6 - S5; */
   if ((err = mp_sub(&S6, &S5, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S1 + S15; */
   if ((err = mp_add(&S1, &S15, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S5 = S5 - t; */
   if ((err = mp_sub(&S5, &t, &S5)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S9 = S9 - S1; */
   if ((err = mp_sub(&S9, &S1, &S9)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S11 = S11 - S1; */
   if ((err = mp_sub(&S11, &S1, &S11)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 - S1; */
   if ((err = mp_sub(&S13, &S1, &S13)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S3 = S3 - S15; */
   if ((err = mp_sub(&S3, &S15, &S3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S4 = S4 - S15; */
   if ((err = mp_sub(&S4, &S15, &S4)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S1<<14; */
   if ((err = mp_mul_2d(&S1, 14, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S4 = S4 - t; */
   if ((err = mp_sub(&S4, &t, &S4)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t<<14; */
   if ((err = mp_mul_2d(&t, 14, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 - t; */
   if ((err = mp_sub(&S3, &t, &S3)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S15<<14; */
   if ((err = mp_mul_2d(&S15, 14, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t<<14; */
   if ((err = mp_mul_2d(&t, 14, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S3 + S11; */
   if ((err = mp_add(&S3, &S11, &S3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S4 = S4 + S9; */
   if ((err = mp_add(&S4, &S9, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 + S10; */
   if ((err = mp_add(&S7, &S10, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S8 = S8 + S12; */
   if ((err = mp_add(&S8, &S12, &S8)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S6<<8; */
   if ((err = mp_mul_2d(&S6, 8, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - t; */
   if ((err = mp_sub(&S7, &t, &S7)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t<<7; */
   if ((err = mp_mul_2d(&t, 7, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S8 = S8 - t; */
   if ((err = mp_sub(&S8, &t, &S8)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t<<7; */
   if ((err = mp_mul_2d(&t, 7, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S2 = S2 + t; */
   if ((err = mp_add(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S4 * 53248; */
   if ((err = mp_mul_d(&S4, 53248, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S5 * 320; */
   if ((err = mp_mul_d(&S5, 320, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S4 = S4 - t; */
   if ((err = mp_sub(&S4, &t, &S4)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S5 * 69632; */
   if ((err = mp_mul_d(&S5, 69632, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S3 = S3 - t; */
   if ((err = mp_sub(&S3, &t, &S3)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S4 * 1360; */
   if ((err = mp_mul_d(&S4, 1360, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S3 = S3 - t; */
   if ((err = mp_sub(&S3, &t, &S3)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = S3 / 11566800; */
   if ((err = mp_div_d(&S3, 11566800, &S3, NULL)) != MP_OKAY)                                            goto LTM_ERR;
   /** t = S3 * 3780; */
   if ((err = mp_mul_d(&S3, 3780, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S4 = S4 - t; */
   if ((err = mp_sub(&S4, &t, &S4)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 / 720; */
   if ((err = mp_div_d(&S4, 720, &S4, NULL)) != MP_OKAY)                                                 goto LTM_ERR;
   /** S5 = S5 - S3; */
   if ((err = mp_sub(&S5, &S3, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 - S4; */
   if ((err = mp_sub(&S5, &S4, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S4 * 1018368000; */
   if ((err = mp_mul_d(&S4, 1018368000, &t)) != MP_OKAY)                                                 goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S3 * 68501160000; */
   if ((err = mp_mul_d(&S3, 68501160000, &t)) != MP_OKAY)                                                goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S2 = -S2; */
   S2.sign = (S2.sign == MP_NEG)?MP_ZPOS: MP_NEG;
   /** t = S7 * 800; */
   if ((err = mp_mul_d(&S7, 800, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S8 = S8 - t; */
   if ((err = mp_sub(&S8, &t, &S8)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S7 * 451584; */
   if ((err = mp_mul_d(&S7, 451584, &t)) != MP_OKAY)                                                     goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S8 * 2856; */
   if ((err = mp_mul_d(&S8, 2856, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S2 = S2 / 372734346600; */
   if ((err = mp_div_d(&S2, 372734346600, &S2, NULL)) != MP_OKAY)                                        goto LTM_ERR;
   /** t = S2 * 7938; */
   if ((err = mp_mul_d(&S2, 7938, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S7 = S7 - t; */
   if ((err = mp_sub(&S7, &t, &S7)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t * 7650; */
   if ((err = mp_mul_d(&t, 7650, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S8 = S8 - t; */
   if ((err = mp_sub(&S8, &t, &S8)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S8 = S8 / 2721600; */
   if ((err = mp_div_d(&S8, 2721600, &S8, NULL)) != MP_OKAY)                                             goto LTM_ERR;
   /** t = S8 * 1800; */
   if ((err = mp_mul_d(&S8, 1800, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S7 = S7 - t; */
   if ((err = mp_sub(&S7, &t, &S7)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S7 = S7 / 288; */
   if ((err = mp_div_d(&S7, 288, &S7, NULL)) != MP_OKAY)                                                 goto LTM_ERR;
   /** S6 = S6 - S2; */
   if ((err = mp_sub(&S6, &S2, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - S8; */
   if ((err = mp_sub(&S6, &S8, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - S7; */
   if ((err = mp_sub(&S6, &S7, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S2<<1; */
   if ((err = mp_mul_2d(&S2, 1, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<1; */
   if ((err = mp_mul_2d(&t, 1, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<1; */
   if ((err = mp_mul_2d(&t, 1, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S8<<3; */
   if ((err = mp_mul_2d(&S8, 3, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<3; */
   if ((err = mp_mul_2d(&t, 3, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<3; */
   if ((err = mp_mul_2d(&t, 3, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S7<<5; */
   if ((err = mp_mul_2d(&S7, 5, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<5; */
   if ((err = mp_mul_2d(&t, 5, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<5; */
   if ((err = mp_mul_2d(&t, 5, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S6<<7; */
   if ((err = mp_mul_2d(&S6, 7, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<7; */
   if ((err = mp_mul_2d(&t, 7, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<7; */
   if ((err = mp_mul_2d(&t, 7, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S3<<2; */
   if ((err = mp_mul_2d(&S3, 2, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t<<2; */
   if ((err = mp_mul_2d(&t, 2, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<2; */
   if ((err = mp_mul_2d(&t, 2, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S4<<4; */
   if ((err = mp_mul_2d(&S4, 4, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t<<4; */
   if ((err = mp_mul_2d(&t, 4, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<4; */
   if ((err = mp_mul_2d(&t, 4, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S5<<6; */
   if ((err = mp_mul_2d(&S5, 6, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t<<6; */
   if ((err = mp_mul_2d(&t, 6, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<6; */
   if ((err = mp_mul_2d(&t, 6, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S9 * 320; */
   if ((err = mp_mul_d(&S9, 320, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S9 * 86016; */
   if ((err = mp_mul_d(&S9, 86016, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S11 * 1360; */
   if ((err = mp_mul_d(&S11, 1360, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S13 = S13 / 47331345600; */
   if ((err = mp_div_d(&S13, 47331345600, &S13, NULL)) != MP_OKAY)                                       goto LTM_ERR;
   /** t = S13 * 15467760; */
   if ((err = mp_mul_d(&S13, 15467760, &t)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S11 = S11 / 725760; */
   if ((err = mp_div_d(&S11, 725760, &S11, NULL)) != MP_OKAY)                                            goto LTM_ERR;
   /** t = S13 * 4092; */
   if ((err = mp_mul_d(&S13, 4092, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S11 * 1008; */
   if ((err = mp_mul_d(&S11, 1008, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S9 = S9 / 192; */
   if ((err = mp_div_d(&S9, 192, &S9, NULL)) != MP_OKAY)                                                 goto LTM_ERR;
   /** t = S10 * 279552; */
   if ((err = mp_mul_d(&S10, 279552, &t)) != MP_OKAY)                                                    goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S10 * 544; */
   if ((err = mp_mul_d(&S10, 544, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S12 * 2600; */
   if ((err = mp_mul_d(&S12, 2600, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S14 = S14 / 384567183000; */
   if ((err = mp_div_d(&S14, 384567183000, &S14, NULL)) != MP_OKAY)                                      goto LTM_ERR;
   /** t = S14 * 62653500; */
   if ((err = mp_mul_d(&S14, 62653500, &t)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S12 = S12 / 3084480; */
   if ((err = mp_div_d(&S12, 3084480, &S12, NULL)) != MP_OKAY)                                           goto LTM_ERR;
   /** t = S14 * 8190; */
   if ((err = mp_mul_d(&S14, 8190, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S12 * 2040; */
   if ((err = mp_mul_d(&S12, 2040, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S10 = S10 / 480; */
   if ((err = mp_div_d(&S10, 480, &S10, NULL)) != MP_OKAY)                                               goto LTM_ERR;
   /** S5 = S5 - S9; */
   if ((err = mp_sub(&S5, &S9, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - S11; */
   if ((err = mp_sub(&S4, &S11, &S4)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S3 - S13; */
   if ((err = mp_sub(&S3, &S13, &S3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S2 = S2 - S14; */
   if ((err = mp_sub(&S2, &S14, &S2)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = S7 - S10; */
   if ((err = mp_sub(&S7, &S10, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S8 = S8 - S12; */
   if ((err = mp_sub(&S8, &S12, &S8)) != MP_OKAY)                                                        goto LTM_ERR;
   /* Needs PARI / GP > 2.6 */
   /** [S4,S8] = [S8,S4]; */
   mp_exch(&S4, &S8);
   /** [S5,S8] = [S8,S5]; */
   mp_exch(&S5, &S8);
   /** [S6,S7] = [S7,S6]; */
   mp_exch(&S6, &S7);
   /** [S7,S8] = [S8,S7]; */
   mp_exch(&S7, &S8);


   /* P  = S1 + S2 * x + S3 * x^2 + S4 * x^3 + S5 * x^4 + S6 * x^5 + S7 * x^6 + S8 * x^7;
      P += S9 * x^8 + S10 * x^9 + S11 * x^10 + S12 * x^11 + S13 * x^12 + S14 * x^13 + S15 * x^14; */

   if ((err = mp_lshd(&S2, 1 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S2, &S1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S3, 2 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S3, &S1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S4, 3 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S4, &S1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S5, 4 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S5, &S1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S6, 5 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S6, &S1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S7, 6 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S7, &S1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S8, 7 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S8, &S1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S9, 8 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S9, &S1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S10, 9 * B)) != MP_OKAY)                                                          goto LTM_ERR;
   if ((err = mp_add(&S10, &S1, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_lshd(&S11, 10 * B)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S11, &S1, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_lshd(&S12, 11 * B)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S12, &S1, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_lshd(&S13, 12 * B)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S13, &S1, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_lshd(&S14, 13 * B)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S14, &S1, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_lshd(&S15, 14 * B)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S15, &S1, c)) != MP_OKAY)                                                          goto LTM_ERR;

   /* P - A^2 */

LTM_ERR:
   mp_clear(&a7);
LTM_ERRa7:
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
   mp_clear_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, &S8, &S9, &S10, &S11, &S12, &S13, &S14, &S15, &t, NULL);
   return err;
}

#endif
