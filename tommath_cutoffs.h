/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
/*
   Current values evaluated on an AMD A8-6600K (64-bit).
   Type "make tune" to optimize them for your machine but
   be aware that it may take a long time. It took 2:30 minutes
   on the aforementioned machine for example.
 */
#define MP_DEFAULT_MUL_KARATSUBA_CUTOFF 105
#define MP_DEFAULT_SQR_KARATSUBA_CUTOFF 154
#define MP_DEFAULT_MUL_TOOM_CUTOFF      131
#define MP_DEFAULT_SQR_TOOM_CUTOFF      192
#define MP_DEFAULT_MUL_TOOM_4_CUTOFF    196
#define MP_DEFAULT_SQR_TOOM_4_CUTOFF    256
#define MP_DEFAULT_MUL_TOOM_5_CUTOFF    239
#define MP_DEFAULT_SQR_TOOM_5_CUTOFF    256
