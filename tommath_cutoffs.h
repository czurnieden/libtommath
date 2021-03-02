/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
/*
   Current values evaluated on an AMD A8-6600K (64-bit).
   Type "make tune" to optimize them for your machine but
   be aware that it may take a long time. It took 2:30 minutes
   on the aforementioned machine for example.
 */
#define MP_DEFAULT_MUL_KARATSUBA_CUTOFF 115
#define MP_DEFAULT_SQR_KARATSUBA_CUTOFF 80
#define MP_DEFAULT_MUL_TOOM_CUTOFF      192
#define MP_DEFAULT_SQR_TOOM_CUTOFF      65
#define MP_DEFAULT_MUL_TOOM_4_CUTOFF    855
#define MP_DEFAULT_SQR_TOOM_4_CUTOFF    108
#define MP_DEFAULT_MUL_TOOM_5_CUTOFF    1061
#define MP_DEFAULT_SQR_TOOM_5_CUTOFF    234
#define MP_DEFAULT_MUL_TOOM_6_CUTOFF    3058
#define MP_DEFAULT_SQR_TOOM_6_CUTOFF    9980
#define MP_DEFAULT_MUL_TOOM_7_CUTOFF    1403
#define MP_DEFAULT_SQR_TOOM_7_CUTOFF    4206
#define MP_DEFAULT_MUL_TOOM_8_CUTOFF    1024
#define MP_DEFAULT_SQR_TOOM_8_CUTOFF    2223
#define MP_DEFAULT_MUL_TOOM_9_CUTOFF    2080
#define MP_DEFAULT_SQR_TOOM_9_CUTOFF    2490
