#include "tommath_private.h"
#ifdef S_MP_SQR_TOOM_5_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */




/*
   This file contains code from J. Arndt's book  "Matters Computational"
   and the accompanying FXT-library with permission of the author.
*/

/*
     Bodrato, Marco, and Alberto Zanoni. "What about Toom-Cook matrices optimality."
     Centro Vito Volterra Universita di Roma Tor Vergata (2006).
*/
mp_err s_mp_sqr_toom_5(const mp_int *a, mp_int *c)
{
   mp_int S1, S2, S3, S4, S5, S6, S7, S8, S9;
   mp_int a0, a1, a2, a3, a4;
   mp_err err = MP_OKAY;
   int B;

   B = a->used / 5;

   if ((err = mp_init_multi(&S1, &S2, &S3, &S4, &S5, &S6, &S7, &S8, &S9, NULL)) != MP_OKAY) {
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


   /** S1 = a0^2 ; */
   if ((err = mp_sqr(&a0, &S1)) != MP_OKAY)                                                              goto LTM_ERR;
   /** S2 = a4^2 ; */
   if ((err = mp_sqr(&a4, &S2)) != MP_OKAY)                                                              goto LTM_ERR;

   /** \\S3 = (a0 + a1 + a2 + a3 + a4)^2; */
   /** S3 = a0 + a1; */
   if ((err = mp_add(&a0, &a1, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 + a2; */
   if ((err = mp_add(&S3, &a2, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 + a3; */
   if ((err = mp_add(&S3, &a3, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 + a4; */
   if ((err = mp_add(&S3, &a4, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3^2; */
   if ((err = mp_sqr(&S3, &S3)) != MP_OKAY)                                                              goto LTM_ERR;

   /** \\S4 = (a0 - a1 + a2 - a3 + a4)^2 ; */
   /** S4 = a0 - a1; */
   if ((err = mp_sub(&a0, &a1, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 + a2; */
   if ((err = mp_add(&S4, &a2, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - a3; */
   if ((err = mp_sub(&S4, &a3, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 + a4; */
   if ((err = mp_add(&S4, &a4, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4^2; */
   if ((err = mp_sqr(&S4, &S4)) != MP_OKAY)                                                              goto LTM_ERR;

   /** \\S5 = 2 * (a0 - a2 + a4) * (a1 - a3); */
   /** \\   = ((a0 - a2 + a4) * (a1 - a3)) << 1; */
   /** S5 = a0 - a2; */
   if ((err = mp_sub(&a0, &a2, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 + a4; */
   if ((err = mp_add(&S5, &a4, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = a1 - a3; */
   if ((err = mp_sub(&a1, &a3, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 * S6; */
   if ((err = mp_mul(&S5, &S6, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 * 2; */
   if ((err = mp_mul_2(&S5, &S5)) != MP_OKAY)                                                            goto LTM_ERR;

   /** \\S6 = (a0 + a1 - a2 - a3 + a4) * (a0 - a1 - a2 + a3 + a4) ; */
   /** S6 = a0 + a1; */
   if ((err = mp_add(&a0, &a1, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - a2; */
   if ((err = mp_sub(&S6, &a2, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - a3; */
   if ((err = mp_sub(&S6, &a3, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 + a4; */
   if ((err = mp_add(&S6, &a4, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = a0 - a1; */
   if ((err = mp_sub(&a0, &a1, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - a2; */
   if ((err = mp_sub(&S7, &a2, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 + a3; */
   if ((err = mp_add(&S7, &a3, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 + a4; */
   if ((err = mp_add(&S7, &a4, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 * S7; */
   if ((err = mp_mul(&S6, &S7, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S7 = (a1 + a2 - a4) * (a1 - a2 - a4 + 2*(a0 - a3) ) ; */
   /** \\S7 = (a1 + a2 - a4) * (2*(a0 - a3) + a1 - a2 - a4  ); */
   /** \\S7 =      S7        *             S8; */
   /** S8 =  a0 - a3; */
   if ((err = mp_sub(&a0, &a3, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = S8 * 2; */
   if ((err = mp_mul_2(&S8, &S8)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S8 = S8 + a1; */
   if ((err = mp_add(&S8, &a1, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = S8 - a2; */
   if ((err = mp_sub(&S8, &a2, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = S8 - a4; */
   if ((err = mp_sub(&S8, &a4, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = a1 + a2; */
   if ((err = mp_add(&a1, &a2, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - a4; */
   if ((err = mp_sub(&S7, &a4, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 * S8; */
   if ((err = mp_mul(&S7, &S8, &S7)) != MP_OKAY)                                                         goto LTM_ERR;

   /** \\S8 = 2 * a0 * a1 ; */
   /** S8 = a0 * a1 ; */
   if ((err = mp_mul(&a0, &a1, &S8)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S8 = S8 * 2; */
   if ((err = mp_mul_2(&S8, &S8)) != MP_OKAY)                                                            goto LTM_ERR;
   /** \\S9 = 2 * a3 * a4 ; */
   /** S9 =  a3 * a4; */
   if ((err = mp_mul(&a3, &a4, &S9)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S9 = S9 * 2; */
   if ((err = mp_mul_2(&S9, &S9)) != MP_OKAY)                                                            goto LTM_ERR;

   /** \\S4 = (S4 + S3) / 2 ; */
   /** S4 = S4 + S3; */
   if ((err = mp_add(&S4, &S3, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 / 2; */
   if ((err = mp_div_2(&S4, &S4)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S3 = S3 - S4 ; */
   if ((err = mp_sub(&S3, &S4, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S6 = (S6 + S4) / 2 ; */
   /** S6 = S6 + S4; */
   if ((err = mp_add(&S6, &S4, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6  / 2 ; */
   if ((err =  mp_div_2(&S6, &S6)) != MP_OKAY)                                                           goto LTM_ERR;
   /** \\S5 = (S3 - S5) / 2 ; */
   /** S5 = S3 - S5; */
   if ((err = mp_sub(&S3, &S5, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 / 2; */
   if ((err = mp_div_2(&S5, &S5)) != MP_OKAY)                                                            goto LTM_ERR;
   /** S4 = S4 - S6 ; */
   if ((err = mp_sub(&S4, &S6, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S3 = S3 - S5 - S8 ; */
   /** S3 = S3 - S5 ; */
   if ((err = mp_sub(&S3, &S5, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S3 = S3 - S8; */
   if ((err = mp_sub(&S3, &S8, &S3)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S6 = S6 - S2 - S1 ; */
   /** S6 = S6 - S2; */
   if ((err = mp_sub(&S6, &S2, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S6 = S6 - S1; */
   if ((err = mp_sub(&S6, &S1, &S6)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S5 = S5 - S9 ; */
   if ((err = mp_sub(&S5, &S9, &S5)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\S7 = S7 - S2 - S8 - S9 + S6 + S3 ; */
   /** S7 = S7 - S2; */
   if ((err = mp_sub(&S7, &S2, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - S8; */
   if ((err = mp_sub(&S7, &S8, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 - S9; */
   if ((err = mp_sub(&S7, &S9, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 + S6; */
   if ((err = mp_add(&S7, &S6, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S7 = S7 + S3; */
   if ((err = mp_add(&S7, &S3, &S7)) != MP_OKAY)                                                         goto LTM_ERR;
   /** S4 = S4 - S7 ; */
   if ((err = mp_sub(&S4, &S7, &S4)) != MP_OKAY)                                                         goto LTM_ERR;
   /** \\ coefficients not in order!; */
   /** P = S2*x^8+ S9*x^7+ S4*x^6+ S3*x^5+ S6*x^4+ S5*x^3+ S7*x^2+ S8*x+ S1 ; */

   if ((err = mp_lshd(&S2, 8 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_lshd(&S9, 7 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S2, &S9, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S4, 6 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S2, &S4, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S3, 5 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S2, &S3, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S6, 4 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S2, &S6, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S5, 3 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S2, &S5, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S7, 2 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S2, &S7, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_lshd(&S8, 1 * B)) != MP_OKAY)                                                           goto LTM_ERR;
   if ((err = mp_add(&S2, &S8, &S2)) != MP_OKAY)                                                         goto LTM_ERR;
   if ((err = mp_add(&S2, &S1, c)) != MP_OKAY)                                                           goto LTM_ERR;

   /** A^2 - P */


LTM_ERR:
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
