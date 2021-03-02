#include "tommath_private.h"
#ifdef S_MP_SQR_TOOM_9_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

mp_err s_mp_sqr_toom_9(const mp_int *a, mp_int *c)
{
   mp_int S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, t;
   mp_int a0, a1, a2, a3, a4, a5, a6, a7, a8;
   mp_err err = MP_OKAY;
   int B;

   B = a->used / 9;

   if ((err = mp_init_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, &S8, &S9, &S10, &S11, &S12, &S13, &S14, &S15, &S16,
                            &S17, &t, NULL)) != MP_OKAY) {
      return err;
   }

   if ((err = mp_init_size(&a0, B)) != MP_OKAY)                                                          goto LTM_ERRa0;
   if ((err = mp_init_size(&a1, B)) != MP_OKAY)                                                          goto LTM_ERRa1;
   if ((err = mp_init_size(&a2, B)) != MP_OKAY)                                                          goto LTM_ERRa2;
   if ((err = mp_init_size(&a3, B)) != MP_OKAY)                                                          goto LTM_ERRa3;
   if ((err = mp_init_size(&a4, B)) != MP_OKAY)                                                          goto LTM_ERRa4;
   if ((err = mp_init_size(&a5, B)) != MP_OKAY)                                                          goto LTM_ERRa5;
   if ((err = mp_init_size(&a6, B)) != MP_OKAY)                                                          goto LTM_ERRa6;
   if ((err = mp_init_size(&a7, B)) != MP_OKAY)                                                          goto LTM_ERRa7;
   if ((err = mp_init_size(&a8, a->used - 8 * B)) != MP_OKAY)                                            goto LTM_ERRa8;

   /* A = a0 * x^0 + a1 * x^1 + a2 * x^2 + a3 * x^3  + a4 * x^4 + a5 * x^5 + a6 * x^6 + a7 * x^7 + a8 * x^8; */
   a0.used = a1.used = a2.used = a3.used = a4.used = a5.used = a6.used = a7.used = B;
   a8.used = a->used - 8 * B;
   s_mp_copy_digs(a0.dp, a->dp, a0.used);
   s_mp_copy_digs(a1.dp, a->dp + B, a1.used);
   s_mp_copy_digs(a2.dp, a->dp + 2 * B, a2.used);
   s_mp_copy_digs(a3.dp, a->dp + 3 * B, a3.used);
   s_mp_copy_digs(a4.dp, a->dp + 4 * B, a4.used);
   s_mp_copy_digs(a5.dp, a->dp + 5 * B, a5.used);
   s_mp_copy_digs(a6.dp, a->dp + 6 * B, a6.used);
   s_mp_copy_digs(a7.dp, a->dp + 7 * B, a7.used);
   s_mp_copy_digs(a8.dp, a->dp + 8 * B, a8.used);
   mp_clamp(&a0);
   mp_clamp(&a1);
   mp_clamp(&a2);
   mp_clamp(&a3);
   mp_clamp(&a4);
   mp_clamp(&a5);
   mp_clamp(&a6);
   mp_clamp(&a7);
   mp_clamp(&a8);




   /** S1 = a0 + a2; */
   if ((err = mp_add(&a0, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + a4; */
   if ((err = mp_add(&S1, &a4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + a6; */
   if ((err = mp_add(&S1, &a6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1 + a8; */
   if ((err = mp_add(&S1, &a8, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S17 = a1 + a3; */
   if ((err = mp_add(&a1, &a3, &S17)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S17 = S17 + a5; */
   if ((err = mp_add(&S17, &a5, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17 + a7; */
   if ((err = mp_add(&S17, &a7, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S2 = S1 - S17; */
   if ((err = mp_sub(&S1, &S17, &S2)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S1 + S17; */
   if ((err = mp_add(&S1, &S17, &S3)) != MP_OKAY)                                                        goto LTM_ERR;



   /** S6 = S2^2; */
   if ((err = mp_sqr(&S2, &S6)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S7 = S3^2; */
   if ((err = mp_sqr(&S3, &S7)) != MP_OKAY)                                                              goto LTM_ERR;



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
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a8; */
   if ((err = mp_add(&S1, &a8, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S17 = a1<<4; */
   if ((err = mp_mul_2d(&a1, 4, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17 + a3; */
   if ((err = mp_add(&S17, &a3, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<4; */
   if ((err = mp_mul_2d(&S17, 4, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a5; */
   if ((err = mp_add(&S17, &a5, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<4; */
   if ((err = mp_mul_2d(&S17, 4, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a7; */
   if ((err = mp_add(&S17, &a7, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<2; */
   if ((err = mp_mul_2d(&S17, 2, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S2 = S1 - S17; */
   if ((err = mp_sub(&S1, &S17, &S2)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S1 + S17; */
   if ((err = mp_add(&S1, &S17, &S3)) != MP_OKAY)                                                        goto LTM_ERR;


   /** S4 = S2^2; */
   if ((err = mp_sqr(&S2, &S4)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S9 = S3^2; */
   if ((err = mp_sqr(&S3, &S9)) != MP_OKAY)                                                              goto LTM_ERR;



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
   /** S1 = S1<<2; */
   if ((err = mp_mul_2d(&S1, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a8; */
   if ((err = mp_add(&S1, &a8, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S17 = a1<<2; */
   if ((err = mp_mul_2d(&a1, 2, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17 + a3; */
   if ((err = mp_add(&S17, &a3, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<2; */
   if ((err = mp_mul_2d(&S17, 2, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a5; */
   if ((err = mp_add(&S17, &a5, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<2; */
   if ((err = mp_mul_2d(&S17, 2, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a7; */
   if ((err = mp_add(&S17, &a7, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<1; */
   if ((err = mp_mul_2d(&S17, 1, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S11 = S1 - S17; */
   if ((err = mp_sub(&S1, &S17, &S11)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S12 = S1 + S17; */
   if ((err = mp_add(&S1, &S17, &S12)) != MP_OKAY)                                                       goto LTM_ERR;


   /** S5 = S11^2; */
   if ((err = mp_sqr(&S11, &S5)) != MP_OKAY)                                                             goto LTM_ERR;
   /** S8 = S12^2; */
   if ((err = mp_sqr(&S12, &S8)) != MP_OKAY)                                                             goto LTM_ERR;


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
   /** S1 = S1<<6; */
   if ((err = mp_mul_2d(&S1, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a8; */
   if ((err = mp_add(&S1, &a8, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S17 = a1<<6; */
   if ((err = mp_mul_2d(&a1, 6, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17 + a3; */
   if ((err = mp_add(&S17, &a3, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<6; */
   if ((err = mp_mul_2d(&S17, 6, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a5; */
   if ((err = mp_add(&S17, &a5, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<6; */
   if ((err = mp_mul_2d(&S17, 6, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a7; */
   if ((err = mp_add(&S17, &a7, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<3; */
   if ((err = mp_mul_2d(&S17, 3, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S12 = S1 - S17; */
   if ((err = mp_sub(&S1, &S17, &S12)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S1 + S17; */
   if ((err = mp_add(&S1, &S17, &S13)) != MP_OKAY)                                                       goto LTM_ERR;

   /** S3 = S12^2; */
   if ((err = mp_sqr(&S12, &S3)) != MP_OKAY)                                                             goto LTM_ERR;
   /** S2 = S13^2; */
   if ((err = mp_sqr(&S13, &S2)) != MP_OKAY)                                                             goto LTM_ERR;



   /** S1 = a8<<2; */
   if ((err = mp_mul_2d(&a8, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a6; */
   if ((err = mp_add(&S1, &a6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<2; */
   if ((err = mp_mul_2d(&S1, 2, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
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
   /** S17 = a7<<2; */
   if ((err = mp_mul_2d(&a7, 2, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17 + a5; */
   if ((err = mp_add(&S17, &a5, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<2; */
   if ((err = mp_mul_2d(&S17, 2, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a3; */
   if ((err = mp_add(&S17, &a3, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<2; */
   if ((err = mp_mul_2d(&S17, 2, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a1; */
   if ((err = mp_add(&S17, &a1, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<1; */
   if ((err = mp_mul_2d(&S17, 1, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S12 = S1 - S17; */
   if ((err = mp_sub(&S1, &S17, &S12)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S1 + S17; */
   if ((err = mp_add(&S1, &S17, &S13)) != MP_OKAY)                                                       goto LTM_ERR;


   /** S10 = S12^2; */
   if ((err = mp_sqr(&S12, &S10)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S11 = S13^2; */
   if ((err = mp_sqr(&S13, &S11)) != MP_OKAY)                                                            goto LTM_ERR;


   /** S1 = a8<<4; */
   if ((err = mp_mul_2d(&a8, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a6; */
   if ((err = mp_add(&S1, &a6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
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
   /** S17 = a7<<4; */
   if ((err = mp_mul_2d(&a7, 4, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17 + a5; */
   if ((err = mp_add(&S17, &a5, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<4; */
   if ((err = mp_mul_2d(&S17, 4, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a3; */
   if ((err = mp_add(&S17, &a3, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<4; */
   if ((err = mp_mul_2d(&S17, 4, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a1; */
   if ((err = mp_add(&S17, &a1, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<2; */
   if ((err = mp_mul_2d(&S17, 2, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S14 = S1 - S17; */
   if ((err = mp_sub(&S1, &S17, &S14)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S15 = S1 + S17; */
   if ((err = mp_add(&S1, &S17, &S15)) != MP_OKAY)                                                       goto LTM_ERR;


   /** S12 = S14^2; */
   if ((err = mp_sqr(&S14, &S12)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S13 = S15^2; */
   if ((err = mp_sqr(&S15, &S13)) != MP_OKAY)                                                            goto LTM_ERR;



   /** S1 = a8<<6; */
   if ((err = mp_mul_2d(&a8, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a6; */
   if ((err = mp_add(&S1, &a6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<6; */
   if ((err = mp_mul_2d(&S1, 6, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
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
   /** S17 = a7<<6; */
   if ((err = mp_mul_2d(&a7, 6, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17 + a5; */
   if ((err = mp_add(&S17, &a5, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<6; */
   if ((err = mp_mul_2d(&S17, 6, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a3; */
   if ((err = mp_add(&S17, &a3, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<6; */
   if ((err = mp_mul_2d(&S17, 6, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S17 = S17 + a1; */
   if ((err = mp_add(&S17, &a1, &S17)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S17 = S17<<3; */
   if ((err = mp_mul_2d(&S17, 3, &S17)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S1 - S17; */
   if ((err = mp_sub(&S1, &S17, &S15)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S16 = S1 + S17; */
   if ((err = mp_add(&S1, &S17, &S16)) != MP_OKAY)                                                       goto LTM_ERR;


   /** S14 = S15^2; */
   if ((err = mp_sqr(&S15, &S14)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S15 = S16^2; */
   if ((err = mp_sqr(&S16, &S15)) != MP_OKAY)                                                            goto LTM_ERR;




   /** S1 = a8<<4; */
   if ((err = mp_mul_2d(&a8, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a7; */
   if ((err = mp_add(&S1, &a7, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a6; */
   if ((err = mp_add(&S1, &a6, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a5; */
   if ((err = mp_add(&S1, &a5, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a4; */
   if ((err = mp_add(&S1, &a4, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a3; */
   if ((err = mp_add(&S1, &a3, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a2; */
   if ((err = mp_add(&S1, &a2, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a1; */
   if ((err = mp_add(&S1, &a1, &S1)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S1 = S1<<4; */
   if ((err = mp_mul_2d(&S1, 4, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S1 = S1 + a0; */
   if ((err = mp_add(&S1, &a0, &S1)) != MP_OKAY)                                                         goto LTM_ERR;


   /** S16 = S1^2; */
   if ((err = mp_sqr(&S1, &S16)) != MP_OKAY)                                                             goto LTM_ERR;

   /** S1 = a0^2; */
   if ((err = mp_sqr(&a0, &S1)) != MP_OKAY)                                                              goto LTM_ERR;

   /** S17 = a8^2; */
   if ((err = mp_sqr(&a8, &S17)) != MP_OKAY)                                                             goto LTM_ERR;







   /* TODO:  interpolation is the same for 9-way multiplication and 9-way square */

   /** S3 = S3 + S2; */
   if ((err = mp_add(&S3, &S2, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3>>1; */
   if ((err = mp_div_2d(&S3, 1, &S3, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S2 = S2 - S3; */
   if ((err = mp_sub(&S2, &S3, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 + S9; */
   if ((err = mp_add(&S4, &S9, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4>>1; */
   if ((err = mp_div_2d(&S4, 1, &S4, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S9 = S9 - S4; */
   if ((err = mp_sub(&S9, &S4, &S9)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 + S8; */
   if ((err = mp_add(&S5, &S8, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5>>1; */
   if ((err = mp_div_2d(&S5, 1, &S5, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S8 = S8 - S5; */
   if ((err = mp_sub(&S8, &S5, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 + S11; */
   if ((err = mp_add(&S10, &S11, &S10)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S10 = S10>>1; */
   if ((err = mp_div_2d(&S10, 1, &S10, NULL)) != MP_OKAY)                                                goto LTM_ERR;
   /** S11 = S11 - S10; */
   if ((err = mp_sub(&S11, &S10, &S11)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S12 = S12 + S13; */
   if ((err = mp_add(&S12, &S13, &S12)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S12 = S12>>1; */
   if ((err = mp_div_2d(&S12, 1, &S12, NULL)) != MP_OKAY)                                                goto LTM_ERR;
   /** S13 = S13 - S12; */
   if ((err = mp_sub(&S13, &S12, &S13)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S14 = S14 + S15; */
   if ((err = mp_add(&S14, &S15, &S14)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S14 = S14>>1; */
   if ((err = mp_div_2d(&S14, 1, &S14, NULL)) != MP_OKAY)                                                goto LTM_ERR;
   /** S15 = S15 - S14; */
   if ((err = mp_sub(&S15, &S14, &S15)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S6 = S6 + S7; */
   if ((err = mp_add(&S6, &S7, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6>>1; */
   if ((err = mp_div_2d(&S6, 1, &S6, NULL)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S7 = S7 - S6; */
   if ((err = mp_sub(&S7, &S6, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S1<<16; */
   if ((err = mp_mul_2d(&S1, 16, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S5 = S5 - t; */
   if ((err = mp_sub(&S5, &t, &S5)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t<<16; */
   if ((err = mp_mul_2d(&t, 16, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - t; */
   if ((err = mp_sub(&S4, &t, &S4)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t<<16; */
   if ((err = mp_mul_2d(&t, 16, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 - t; */
   if ((err = mp_sub(&S3, &t, &S3)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S6 = S6 - S1; */
   if ((err = mp_sub(&S6, &S1, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 - S1; */
   if ((err = mp_sub(&S10, &S1, &S10)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S12 = S12 - S1; */
   if ((err = mp_sub(&S12, &S1, &S12)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S14 = S14 - S1; */
   if ((err = mp_sub(&S14, &S1, &S14)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S16 = S16 - S1; */
   if ((err = mp_sub(&S16, &S1, &S16)) != MP_OKAY)                                                       goto LTM_ERR;
   /** t = S17<<16; */
   if ((err = mp_mul_2d(&S17, 16, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<16; */
   if ((err = mp_mul_2d(&t, 16, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<16; */
   if ((err = mp_mul_2d(&t, 16, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<16; */
   if ((err = mp_mul_2d(&t, 16, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S3 - S17; */
   if ((err = mp_sub(&S3, &S17, &S3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S4 = S4 - S17; */
   if ((err = mp_sub(&S4, &S17, &S4)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S5 = S5 - S17; */
   if ((err = mp_sub(&S5, &S17, &S5)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S6 = S6 - S17; */
   if ((err = mp_sub(&S6, &S17, &S6)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S2 = S2 + S15; */
   if ((err = mp_add(&S2, &S15, &S2)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S3 + S14; */
   if ((err = mp_add(&S3, &S14, &S3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S4 = S4 + S12; */
   if ((err = mp_add(&S4, &S12, &S4)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S5 = S5 + S10; */
   if ((err = mp_add(&S5, &S10, &S5)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S8 = S8 + S11; */
   if ((err = mp_add(&S8, &S11, &S8)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S9 = S9 + S13; */
   if ((err = mp_add(&S9, &S13, &S9)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S6<<9; */
   if ((err = mp_mul_2d(&S6, 9, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 - t; */
   if ((err = mp_sub(&S5, &t, &S5)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t<<8; */
   if ((err = mp_mul_2d(&t, 8, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 - t; */
   if ((err = mp_sub(&S4, &t, &S4)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = t<<8; */
   if ((err = mp_mul_2d(&t, 8, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = S3 - t; */
   if ((err = mp_sub(&S3, &t, &S3)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S8 * 212992; */
   if ((err = mp_mul_d(&S8, 212992, &t)) != MP_OKAY)                                                     goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S5 * 277022736; */
   if ((err = mp_mul_d(&S5, 277022736, &t)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S3 = S3 - t; */
   if ((err = mp_sub(&S3, &t, &S3)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = -S3; */
   S3.sign = (S3.sign == MP_NEG)?MP_ZPOS: MP_NEG;
   /** t = S5 * 16900; */
   if ((err = mp_mul_d(&S5, 16900, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S4 = S4 - t; */
   if ((err = mp_sub(&S4, &t, &S4)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = -S4; */
   S4.sign = (S4.sign == MP_NEG)?MP_ZPOS: MP_NEG;
   /** t = S4 * 17988; */
   if ((err = mp_mul_d(&S4, 17988, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S3 = S3 - t; */
   if ((err = mp_sub(&S3, &t, &S3)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S3 = S3 / 133641446400; */
   if ((err = mp_div_d(&S3, 133641446400, &S3, NULL)) != MP_OKAY)                                        goto LTM_ERR;
   /** t = S3 * 44193600; */
   if ((err = mp_mul_d(&S3, 44193600, &t)) != MP_OKAY)                                                   goto LTM_ERR;
   /** S4 = S4 - t; */
   if ((err = mp_sub(&S4, &t, &S4)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S4 = S4 / 8812800; */
   if ((err = mp_div_d(&S4, 8812800, &S4, NULL)) != MP_OKAY)                                             goto LTM_ERR;
   /** t = S3 * 3600; */
   if ((err = mp_mul_d(&S3, 3600, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S5 = S5 - t; */
   if ((err = mp_sub(&S5, &t, &S5)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S4 * 576; */
   if ((err = mp_mul_d(&S4, 576, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S5 = S5 - t; */
   if ((err = mp_sub(&S5, &t, &S5)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S5 = S5 / 15876; */
   if ((err = mp_div_d(&S5, 15876, &S5, NULL)) != MP_OKAY)                                               goto LTM_ERR;
   /** S6 = S6 - S5; */
   if ((err = mp_sub(&S6, &S5, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - S4; */
   if ((err = mp_sub(&S6, &S4, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - S3; */
   if ((err = mp_sub(&S6, &S3, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S7 * 640; */
   if ((err = mp_mul_d(&S7, 640, &t)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S8 = S8 - t; */
   if ((err = mp_sub(&S8, &t, &S8)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S8 * 5657600; */
   if ((err = mp_mul_d(&S8, 5657600, &t)) != MP_OKAY)                                                    goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S7 * 278528; */
   if ((err = mp_mul_d(&S7, 278528, &t)) != MP_OKAY)                                                     goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S8 * 2720; */
   if ((err = mp_mul_d(&S8, 2720, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S9 * 35490; */
   if ((err = mp_mul_d(&S9, 35490, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S2 = S2 - t; */
   if ((err = mp_sub(&S2, &t, &S2)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S2 = -S2; */
   S2.sign = (S2.sign == MP_NEG)?MP_ZPOS: MP_NEG;
   /** S2 = S2 / 1136785104000; */
   if ((err = mp_div_d(&S2, 1136785104000, &S2, NULL)) != MP_OKAY)                                       goto LTM_ERR;
   /** t = S2 * 46267200; */
   if ((err = mp_mul_d(&S2, 46267200, &t)) != MP_OKAY)                                                   goto LTM_ERR;
   /** S9 = S9 - t; */
   if ((err = mp_sub(&S9, &t, &S9)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S9 = S9 / 986069700; */
   if ((err = mp_div_d(&S9, 986069700, &S9, NULL)) != MP_OKAY)                                           goto LTM_ERR;
   /** t = S9 * 32130; */
   if ((err = mp_mul_d(&S9, 32130, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S8 = S8 - t; */
   if ((err = mp_sub(&S8, &t, &S8)) != MP_OKAY)                                                          goto LTM_ERR;
   /** t = S2 * 7560; */
   if ((err = mp_mul_d(&S2, 7560, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S8 = S8 - t; */
   if ((err = mp_sub(&S8, &t, &S8)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S8 = S8 / 1440; */
   if ((err = mp_div_d(&S8, 1440, &S8, NULL)) != MP_OKAY)                                                goto LTM_ERR;
   /** S7 = S7 - S9; */
   if ((err = mp_sub(&S7, &S9, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - S2; */
   if ((err = mp_sub(&S7, &S2, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - S8; */
   if ((err = mp_sub(&S7, &S8, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** t = S5<<2; */
   if ((err = mp_mul_2d(&S5, 2, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<2; */
   if ((err = mp_mul_2d(&t, 2, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<2; */
   if ((err = mp_mul_2d(&t, 2, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<2; */
   if ((err = mp_mul_2d(&t, 2, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S9<<1; */
   if ((err = mp_mul_2d(&S9, 1, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<1; */
   if ((err = mp_mul_2d(&t, 1, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<1; */
   if ((err = mp_mul_2d(&t, 1, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S15 = S15 - t; */
   if ((err = mp_sub(&S15, &t, &S15)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<1; */
   if ((err = mp_mul_2d(&t, 1, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S2<<3; */
   if ((err = mp_mul_2d(&S2, 3, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<3; */
   if ((err = mp_mul_2d(&t, 3, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<3; */
   if ((err = mp_mul_2d(&t, 3, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S15 = S15 - t; */
   if ((err = mp_sub(&S15, &t, &S15)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<3; */
   if ((err = mp_mul_2d(&t, 3, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S3<<4; */
   if ((err = mp_mul_2d(&S3, 4, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<4; */
   if ((err = mp_mul_2d(&t, 4, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<4; */
   if ((err = mp_mul_2d(&t, 4, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<4; */
   if ((err = mp_mul_2d(&t, 4, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S8<<5; */
   if ((err = mp_mul_2d(&S8, 5, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<5; */
   if ((err = mp_mul_2d(&t, 5, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<5; */
   if ((err = mp_mul_2d(&t, 5, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S15 = S15 - t; */
   if ((err = mp_sub(&S15, &t, &S15)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<5; */
   if ((err = mp_mul_2d(&t, 5, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S4<<6; */
   if ((err = mp_mul_2d(&S4, 6, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<6; */
   if ((err = mp_mul_2d(&t, 6, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<6; */
   if ((err = mp_mul_2d(&t, 6, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<6; */
   if ((err = mp_mul_2d(&t, 6, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S7<<7; */
   if ((err = mp_mul_2d(&S7, 7, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<7; */
   if ((err = mp_mul_2d(&t, 7, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<7; */
   if ((err = mp_mul_2d(&t, 7, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S15 = S15 - t; */
   if ((err = mp_sub(&S15, &t, &S15)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<7; */
   if ((err = mp_mul_2d(&t, 7, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S6<<8; */
   if ((err = mp_mul_2d(&S6, 8, &t)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<8; */
   if ((err = mp_mul_2d(&t, 8, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<8; */
   if ((err = mp_mul_2d(&t, 8, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = t<<8; */
   if ((err = mp_mul_2d(&t, 8, &t)) != MP_OKAY)                                                          goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S13 * 278528; */
   if ((err = mp_mul_d(&S13, 278528, &t)) != MP_OKAY)                                                    goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S11 * 344064; */
   if ((err = mp_mul_d(&S11, 344064, &t)) != MP_OKAY)                                                    goto LTM_ERR;
   /** S15 = S15 - t; */
   if ((err = mp_sub(&S15, &t, &S15)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S11 * 640; */
   if ((err = mp_mul_d(&S11, 640, &t)) != MP_OKAY)                                                       goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S12 * 1052672; */
   if ((err = mp_mul_d(&S12, 1052672, &t)) != MP_OKAY)                                                   goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S10 * 1118208; */
   if ((err = mp_mul_d(&S10, 1118208, &t)) != MP_OKAY)                                                   goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S10 * 1088; */
   if ((err = mp_mul_d(&S10, 1088, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S13 * 5657600; */
   if ((err = mp_mul_d(&S13, 5657600, &t)) != MP_OKAY)                                                   goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S13 * 2720; */
   if ((err = mp_mul_d(&S13, 2720, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S15 = S15 - t; */
   if ((err = mp_sub(&S15, &t, &S15)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S14 * 4112; */
   if ((err = mp_mul_d(&S14, 4112, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S12 * 5200; */
   if ((err = mp_mul_d(&S12, 5200, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S14 = S14 - t; */
   if ((err = mp_sub(&S14, &t, &S14)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S14 * 17476; */
   if ((err = mp_mul_d(&S14, 17476, &t)) != MP_OKAY)                                                     goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S14 = S14 / 3076537464000; */
   if ((err = mp_div_d(&S14, 3076537464000, &S14, NULL)) != MP_OKAY)                                     goto LTM_ERR;
   /** t = S14 * 250614000; */
   if ((err = mp_mul_d(&S14, 250614000, &t)) != MP_OKAY)                                                 goto LTM_ERR;
   /** S12 = S12 - t; */
   if ((err = mp_sub(&S12, &t, &S12)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S12 = S12 / 12337920; */
   if ((err = mp_div_d(&S12, 12337920, &S12, NULL)) != MP_OKAY)                                          goto LTM_ERR;
   /** t = S14 * 16380; */
   if ((err = mp_mul_d(&S14, 16380, &t)) != MP_OKAY)                                                     goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S12 * 4080; */
   if ((err = mp_mul_d(&S12, 4080, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S10 = S10 - t; */
   if ((err = mp_sub(&S10, &t, &S10)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S10 = S10 / 960; */
   if ((err = mp_div_d(&S10, 960, &S10, NULL)) != MP_OKAY)                                               goto LTM_ERR;
   /** t = S15 * 35490; */
   if ((err = mp_mul_d(&S15, 35490, &t)) != MP_OKAY)                                                     goto LTM_ERR;
   /** S16 = S16 - t; */
   if ((err = mp_sub(&S16, &t, &S16)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S16 = -S16; */
   S16.sign = (S16.sign == MP_NEG)?MP_ZPOS: MP_NEG;
   /** S16 = S16 / 9303449291136000; */
   if ((err = mp_div_d(&S16, 9303449291136000, &S16, NULL)) != MP_OKAY)                                  goto LTM_ERR;
   /** t = S16 * 378650764800; */
   if ((err = mp_mul_d(&S16, 378650764800, &t)) != MP_OKAY)                                              goto LTM_ERR;
   /** S15 = S15 - t; */
   if ((err = mp_sub(&S15, &t, &S15)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S15 = S15 / 32309559790200; */
   if ((err = mp_div_d(&S15, 32309559790200, &S15, NULL)) != MP_OKAY)                                    goto LTM_ERR;
   /** t = S15 * 1052771580; */
   if ((err = mp_mul_d(&S15, 1052771580, &t)) != MP_OKAY)                                                goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S15 * 32766; */
   if ((err = mp_mul_d(&S15, 32766, &t)) != MP_OKAY)                                                     goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S16 * 61871040; */
   if ((err = mp_mul_d(&S16, 61871040, &t)) != MP_OKAY)                                                  goto LTM_ERR;
   /** S13 = S13 - t; */
   if ((err = mp_sub(&S13, &t, &S13)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S13 = S13 / 2903040; */
   if ((err = mp_div_d(&S13, 2903040, &S13, NULL)) != MP_OKAY)                                           goto LTM_ERR;
   /** t = S16 * 8184; */
   if ((err = mp_mul_d(&S16, 8184, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** t = S13 * 2016; */
   if ((err = mp_mul_d(&S13, 2016, &t)) != MP_OKAY)                                                      goto LTM_ERR;
   /** S11 = S11 - t; */
   if ((err = mp_sub(&S11, &t, &S11)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S11 = S11 / 384; */
   if ((err = mp_div_d(&S11, 384, &S11, NULL)) != MP_OKAY)                                               goto LTM_ERR;
   /** S2 = S2 - S16; */
   if ((err = mp_sub(&S2, &S16, &S2)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S3 = S3 - S12; */
   if ((err = mp_sub(&S3, &S12, &S3)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S4 = S4 - S10; */
   if ((err = mp_sub(&S4, &S10, &S4)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S5 = S5 - S14; */
   if ((err = mp_sub(&S5, &S14, &S5)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S7 = S7 - S11; */
   if ((err = mp_sub(&S7, &S11, &S7)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S8 = S8 - S13; */
   if ((err = mp_sub(&S8, &S13, &S8)) != MP_OKAY)                                                        goto LTM_ERR;
   /** S9 = S9 - S15; */
   if ((err = mp_sub(&S9, &S15, &S9)) != MP_OKAY)                                                        goto LTM_ERR;
   /** [S2,S9]=[S9,S2] */
   mp_exch(&S2,&S9);
   /** [S3,S5]=[S5,S3] */
   mp_exch(&S3,&S5);
   /** [S4,S9]=[S9,S4] */
   mp_exch(&S4,&S9);
   /** [S6,S8]=[S8,S6] */
   mp_exch(&S6,&S8);
   /** [S7,S9]=[S9,S7] */
   mp_exch(&S7,&S9);
   /** [S8,S9]=[S9,S8] */
   mp_exch(&S8,&S9);
   /** [S10,S11]=[S11,S10] */
   mp_exch(&S10,&S11);
   /** [S12,S13]=[S13,S12] */
   mp_exch(&S12,&S13);
   /** [S14,S16]=[S16,S14] */
   mp_exch(&S14,&S16);
   /** [S15,S16]=[S16,S15] */
   mp_exch(&S15,&S16);

   /** P  = S1 + S2 * x + S3 * x^2 + S4 * x^3 + S5 * x^4 + S6 * x^5 + S7 * x^6 + S8 * x^7;*/
   /** P += S9 * x^8 + S10 * x^9 + S11 * x^10 + S12 * x^11 + S13 * x^12 + S14 * x^13 + S15 * x^14;*/
   /** P += S16 * x^15 + S17 * x^16;*/

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
   if ((err = mp_add(&S15, &S1, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_lshd(&S16, 15 * B)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S16, &S1, &S1)) != MP_OKAY)                                                        goto LTM_ERR;
   if ((err = mp_lshd(&S17, 16 * B)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S17, &S1, c)) != MP_OKAY)                                                          goto LTM_ERR;

   /** P - A*B */

LTM_ERR:
   mp_clear(&a8);
LTM_ERRa8:
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
   mp_clear_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, &S8, &S9, &S10, &S11, &S12, &S13, &S14, &S15, &S16, &S17, &t, NULL);
   return err;
}

#endif
