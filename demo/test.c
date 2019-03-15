#include "shared.h"

static int test_trivial_stuff(void)
{
   mp_int a, b, c, d;
   if (mp_init_multi(&a, &b, &c, &d, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   /* a: 0->5 */
   mp_set_int(&a, 5uL);
   /* a: 5-> b: -5 */
   mp_neg(&a, &b);
   if (mp_cmp(&a, &b) != MP_GT) {
      goto LBL_ERR;
   }
   if (mp_cmp(&b, &a) != MP_LT) {
      goto LBL_ERR;
   }
   /* a: 5-> a: -5 */
   mp_neg(&a, &a);
   if (mp_cmp(&b, &a) != MP_EQ) {
      goto LBL_ERR;
   }
   /* a: -5-> b: 5 */
   mp_abs(&a, &b);
   if (mp_isneg(&b) != MP_NO) {
      goto LBL_ERR;
   }
   /* a: -5-> b: -4 */
   mp_add_d(&a, 1uL, &b);
   if (mp_isneg(&b) != MP_YES) {
      goto LBL_ERR;
   }
   if (mp_get_int(&b) != 4) {
      goto LBL_ERR;
   }
   /* a: -5-> b: 1 */
   mp_add_d(&a, 6uL, &b);
   if (mp_get_int(&b) != 1) {
      goto LBL_ERR;
   }
   /* a: -5-> a: 1 */
   mp_add_d(&a, 6uL, &a);
   if (mp_get_int(&a) != 1) {
      goto LBL_ERR;
   }
   mp_zero(&a);
   /* a: 0-> a: 6 */
   mp_add_d(&a, 6uL, &a);
   if (mp_get_int(&a) != 6) {
      goto LBL_ERR;
   }

   mp_set_int(&a, 42uL);
   mp_set_int(&b, 1uL);
   mp_neg(&b, &b);
   mp_set_int(&c, 1uL);
   mp_exptmod(&a, &b, &c, &d);

   mp_set_int(&c, 7uL);
   mp_exptmod(&a, &b, &c, &d);

   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_FAILURE;
}

static int test_mp_jacobi(void)
{
   struct mp_jacobi_st {
      unsigned long n;
      int c[16];
   };

   static struct mp_jacobi_st jacobi[] = {
      { 3, {  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1 } },
      { 5, {  0,  1, -1, -1,  1,  0,  1, -1, -1,  1,  0,  1, -1, -1,  1,  0 } },
      { 7, {  1, -1,  1, -1, -1,  0,  1,  1, -1,  1, -1, -1,  0,  1,  1, -1 } },
      { 9, { -1,  1,  0,  1,  1,  0,  1,  1,  0,  1,  1,  0,  1,  1,  0,  1 } },
   };

   int i, n, err, should, cnt;

   mp_int a, b;
   if (mp_init_multi(&a, &b, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   mp_set_int(&a, 0uL);
   mp_set_int(&b, 1uL);
   if ((err = mp_jacobi(&a, &b, &i)) != MP_OKAY) {
      printf("Failed executing mp_jacobi(0 | 1) %s.\n", mp_error_to_string(err));
      goto LBL_ERR;
   }
   if (i != 1) {
      printf("Failed trivial mp_jacobi(0 | 1) %d != 1\n", i);
      goto LBL_ERR;
   }
   for (cnt = 0; cnt < (int)(sizeof(jacobi)/sizeof(jacobi[0])); ++cnt) {
      mp_set_int(&b, jacobi[cnt].n);
      /* only test positive values of a */
      for (n = -5; n <= 10; ++n) {
         mp_set_int(&a, abs(n));
         should = MP_OKAY;
         if (n < 0) {
            mp_neg(&a, &a);
            /* Until #44 is fixed the negative a's must fail */
            should = MP_VAL;
         }
         if ((err = mp_jacobi(&a, &b, &i)) != should) {
            printf("Failed executing mp_jacobi(%d | %lu) %s.\n", n, jacobi[cnt].n, mp_error_to_string(err));
            goto LBL_ERR;
         }
         if ((err == MP_OKAY) && (i != jacobi[cnt].c[n + 5])) {
            printf("Failed trivial mp_jacobi(%d | %lu) %d != %d\n", n, jacobi[cnt].n, i, jacobi[cnt].c[n + 5]);
            goto LBL_ERR;
         }
      }
   }

   mp_clear_multi(&a, &b, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, NULL);
   return EXIT_FAILURE;
}

static int test_mp_kronecker(void)
{
   struct mp_kronecker_st {
      long n;
      int c[21];
   };
   static struct mp_kronecker_st kronecker[] = {
      /*-10, -9, -8, -7,-6, -5, -4, -3, -2, -1, 0, 1,  2,  3, 4,  5,  6,  7,  8, 9, 10*/
      { -10, {  0, -1,  0, -1, 0,  0,  0,  1,  0, -1, 0, 1,  0, -1, 0,  0,  0,  1,  0, 1,  0  } },
      {  -9, { -1,  0, -1,  1, 0, -1, -1,  0, -1, -1, 0, 1,  1,  0, 1,  1,  0, -1,  1, 0,  1  } },
      {  -8, {  0, -1,  0,  1, 0,  1,  0, -1,  0, -1, 0, 1,  0,  1, 0, -1,  0, -1,  0, 1,  0  } },
      {  -7, {  1, -1, -1,  0, 1,  1, -1,  1, -1, -1, 0, 1,  1, -1, 1, -1, -1,  0,  1, 1, -1  } },
      {  -6, {  0,  0,  0, -1, 0, -1,  0,  0,  0, -1, 0, 1,  0,  0, 0,  1,  0,  1,  0, 0,  0  } },
      {  -5, {  0, -1,  1, -1, 1,  0, -1, -1,  1, -1, 0, 1, -1,  1, 1,  0, -1,  1, -1, 1,  0  } },
      {  -4, {  0, -1,  0,  1, 0, -1,  0,  1,  0, -1, 0, 1,  0, -1, 0,  1,  0, -1,  0, 1,  0  } },
      {  -3, { -1,  0,  1, -1, 0,  1, -1,  0,  1, -1, 0, 1, -1,  0, 1, -1,  0,  1, -1, 0,  1  } },
      {  -2, {  0, -1,  0,  1, 0,  1,  0, -1,  0, -1, 0, 1,  0,  1, 0, -1,  0, -1,  0, 1,  0  } },
      {  -1, { -1, -1, -1,  1, 1, -1, -1,  1, -1, -1, 1, 1,  1, -1, 1,  1, -1, -1,  1, 1,  1  } },
      {   0, {  0,  0,  0,  0, 0,  0,  0,  0,  0,  1, 0, 1,  0,  0, 0,  0,  0,  0,  0, 0,  0  } },
      {   1, {  1,  1,  1,  1, 1,  1,  1,  1,  1,  1, 1, 1,  1,  1, 1,  1,  1,  1,  1, 1,  1  } },
      {   2, {  0,  1,  0,  1, 0, -1,  0, -1,  0,  1, 0, 1,  0, -1, 0, -1,  0,  1,  0, 1,  0  } },
      {   3, {  1,  0, -1, -1, 0, -1,  1,  0, -1,  1, 0, 1, -1,  0, 1, -1,  0, -1, -1, 0,  1  } },
      {   4, {  0,  1,  0,  1, 0,  1,  0,  1,  0,  1, 0, 1,  0,  1, 0,  1,  0,  1,  0, 1,  0  } },
      {   5, {  0,  1, -1, -1, 1,  0,  1, -1, -1,  1, 0, 1, -1, -1, 1,  0,  1, -1, -1, 1,  0  } },
      {   6, {  0,  0,  0, -1, 0,  1,  0,  0,  0,  1, 0, 1,  0,  0, 0,  1,  0, -1,  0, 0,  0  } },
      {   7, { -1,  1,  1,  0, 1, -1,  1,  1,  1,  1, 0, 1,  1,  1, 1, -1,  1,  0,  1, 1, -1  } },
      {   8, {  0,  1,  0,  1, 0, -1,  0, -1,  0,  1, 0, 1,  0, -1, 0, -1,  0,  1,  0, 1,  0  } },
      {   9, {  1,  0,  1,  1, 0,  1,  1,  0,  1,  1, 0, 1,  1,  0, 1,  1,  0,  1,  1, 0,  1  } },
      {  10, {  0,  1,  0, -1, 0,  0,  0,  1,  0,  1, 0, 1,  0,  1, 0,  0,  0, -1,  0, 1,  0  } }
   };

   long k, m;
   int i, err, cnt;

   mp_int a, b;
   if (mp_init_multi(&a, &b, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   mp_set_int(&a, 0uL);
   mp_set_int(&b, 1uL);
   if ((err = mp_kronecker(&a, &b, &i)) != MP_OKAY) {
      printf("Failed executing mp_kronecker(0 | 1) %s.\n", mp_error_to_string(err));
      goto LBL_ERR;
   }
   if (i != 1) {
      printf("Failed trivial mp_kronecker(0 | 1) %d != 1\n", i);
      goto LBL_ERR;
   }
   for (cnt = 0; cnt < (int)(sizeof(kronecker)/sizeof(kronecker[0])); ++cnt) {
      k = kronecker[cnt].n;
      if (k < 0) {
         mp_set_int(&a, (unsigned long)(-k));
         mp_neg(&a, &a);
      } else {
         mp_set_int(&a, (unsigned long) k);
      }
      /* only test positive values of a */
      for (m = -10; m <= 10; m++) {
         if (m < 0) {
            mp_set_int(&b,(unsigned long)(-m));
            mp_neg(&b, &b);
         } else {
            mp_set_int(&b, (unsigned long) m);
         }
         if ((err = mp_kronecker(&a, &b, &i)) != MP_OKAY) {
            printf("Failed executing mp_kronecker(%ld | %ld) %s.\n", kronecker[cnt].n, m, mp_error_to_string(err));
            goto LBL_ERR;
         }
         if ((err == MP_OKAY) && (i != kronecker[cnt].c[m + 10])) {
            printf("Failed trivial mp_kronecker(%ld | %ld) %d != %d\n", kronecker[cnt].n, m, i, kronecker[cnt].c[m + 10]);
            goto LBL_ERR;
         }
      }
   }

   mp_clear_multi(&a, &b, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, NULL);
   return EXIT_FAILURE;
}

static int test_mp_complement(void)
{
   int i;

   mp_int a, b, c;
   if (mp_init_multi(&a, &b, &c, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   for (i = 0; i < 1000; ++i) {
      int l = (rand() * rand() + 1) * (rand() % 1 ? -1 : 1);
      mp_set_int(&a, labs(l));
      if (l < 0)
         mp_neg(&a, &a);
      mp_complement(&a, &b);

      l = ~l;
      mp_set_int(&c, labs(l));
      if (l < 0)
         mp_neg(&c, &c);

      if (mp_cmp(&b, &c) != MP_EQ) {
         printf("\nmp_complement() bad result!");
         goto LBL_ERR;
      }
   }

   mp_clear_multi(&a, &b, &c, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, NULL);
   return EXIT_FAILURE;
}

static int test_mp_tc_div_2d(void)
{
   int i;

   mp_int a, b, d;
   if (mp_init_multi(&a, &b, &d, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   for (i = 0; i < 1000; ++i) {
      int l, em;

      l = (rand() * rand() + 1) * (rand() % 1 ? -1 : 1);
      mp_set_int(&a, labs(l));
      if (l < 0)
         mp_neg(&a, &a);

      em = rand() % 32;

      mp_set_int(&d, labs(l >> em));
      if ((l >> em) < 0)
         mp_neg(&d, &d);

      mp_tc_div_2d(&a, em, &b);
      if (mp_cmp(&b, &d) != MP_EQ) {
         printf("\nmp_tc_div_2d() bad result!");
         goto LBL_ERR;
      }
   }

   mp_clear_multi(&a, &b, &d, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &d, NULL);
   return EXIT_FAILURE;

}

static int test_mp_tc_xor(void)
{
   int i;

   mp_int a, b, c, d;
   if (mp_init_multi(&a, &b, &c, &d, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   for (i = 0; i < 1000; ++i) {
      int l, em;

      l = (rand() * rand() + 1) * (rand() % 1 ? -1 : 1);
      mp_set_int(&a, labs(l));
      if (l < 0)
         mp_neg(&a, &a);

      em = (rand() * rand() + 1) * (rand() % 1 ? -1 : 1);
      mp_set_int(&b, labs(em));
      if (em < 0)
         mp_neg(&b, &b);

      mp_set_int(&d, labs(l ^ em));
      if ((l ^ em) < 0)
         mp_neg(&d, &d);

      mp_tc_xor(&a, &b, &c);
      if (mp_cmp(&c, &d) != MP_EQ) {
         printf("\nmp_tc_xor() bad result!");
         goto LBL_ERR;
      }
   }

   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_FAILURE;

}

static int test_mp_tc_or(void)
{
   int i;

   mp_int a, b, c, d;
   if (mp_init_multi(&a, &b, &c, &d, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   for (i = 0; i < 1000; ++i) {
      int l, em;

      l = (rand() * rand() + 1) * (rand() % 1 ? -1 : 1);
      mp_set_int(&a, labs(l));
      if (l < 0)
         mp_neg(&a, &a);

      em = (rand() * rand() + 1) * (rand() % 1 ? -1 : 1);
      mp_set_int(&b, labs(em));
      if (em < 0)
         mp_neg(&b, &b);

      mp_set_int(&d, labs(l | em));
      if ((l | em) < 0)
         mp_neg(&d, &d);

      mp_tc_or(&a, &b, &c);
      if (mp_cmp(&c, &d) != MP_EQ) {
         printf("\nmp_tc_or() bad result!");
         goto LBL_ERR;
      }
   }

   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_FAILURE;
}

static int test_mp_tc_and(void)
{
   int i;

   mp_int a, b, c, d;
   if (mp_init_multi(&a, &b, &c, &d, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   for (i = 0; i < 1000; ++i) {
      int l, em;

      l = (rand() * rand() + 1) * (rand() % 1 ? -1 : 1);
      mp_set_int(&a, labs(l));
      if (l < 0)
         mp_neg(&a, &a);

      em = (rand() * rand() + 1) * (rand() % 1 ? -1 : 1);
      mp_set_int(&b, labs(em));
      if (em < 0)
         mp_neg(&b, &b);

      mp_set_int(&d, labs(l & em));
      if ((l & em) < 0)
         mp_neg(&d, &d);

      mp_tc_and(&a, &b, &c);
      if (mp_cmp(&c, &d) != MP_EQ) {
         printf("\nmp_tc_and() bad result!");
         goto LBL_ERR;
      }
   }

   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_FAILURE;
}

static int test_mp_invmod(void)
{
   mp_int a, b, c, d;
   if (mp_init_multi(&a, &b, &c, &d, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   /* mp_invmod corner-case of https://github.com/libtom/libtommath/issues/118 */
   {
      const char *a_ = "47182BB8DF0FFE9F61B1F269BACC066B48BA145D35137D426328DC3F88A5EA44";
      const char *b_ = "FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF";
      const char *should_ = "0521A82E10376F8E4FDEF9A32A427AC2A0FFF686E00290D39E3E4B5522409596";

      if (mp_read_radix(&a, a_, 16) != MP_OKAY) {
         printf("\nmp_read_radix(a) failed!");
         goto LBL_ERR;
      }
      if (mp_read_radix(&b, b_, 16) != MP_OKAY) {
         printf("\nmp_read_radix(b) failed!");
         goto LBL_ERR;
      }
      if (mp_read_radix(&c, should_, 16) != MP_OKAY) {
         printf("\nmp_read_radix(should) failed!");
         goto LBL_ERR;
      }

      if (mp_invmod(&a, &b, &d) != MP_OKAY) {
         printf("\nmp_invmod() failed!");
         goto LBL_ERR;
      }

      if (mp_cmp(&c, &d) != MP_EQ) {
         printf("\nmp_invmod() bad result!");
         goto LBL_ERR;
      }
   }

   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_FAILURE;

}

static int test_mp_set_double(void)
{
   int i;

   mp_int a, b;
   if (mp_init_multi(&a, &b, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   /* test mp_get_double/mp_set_double */
#if defined(__STDC_IEC_559__) || defined(__GCC_IEC_559)
   if (mp_set_double(&a, +1.0/0.0) != MP_VAL) {
      printf("\nmp_set_double should return MP_VAL for +inf");
      goto LBL_ERR;
   }
   if (mp_set_double(&a, -1.0/0.0) != MP_VAL) {
      printf("\nmp_set_double should return MP_VAL for -inf");
      goto LBL_ERR;
   }
   if (mp_set_double(&a, +0.0/0.0) != MP_VAL) {
      printf("\nmp_set_double should return MP_VAL for NaN");
      goto LBL_ERR;
   }
   if (mp_set_double(&a, -0.0/0.0) != MP_VAL) {
      printf("\nmp_set_double should return MP_VAL for NaN");
      goto LBL_ERR;
   }

   for (i = 0; i < 1000; ++i) {
      int tmp = rand();
      double dbl = (double)tmp * rand() + 1;
      if (mp_set_double(&a, dbl) != MP_OKAY) {
         printf("\nmp_set_double() failed");
         goto LBL_ERR;
      }
      if (dbl != mp_get_double(&a)) {
         printf("\nmp_get_double() bad result!");
         goto LBL_ERR;
      }
      if (mp_set_double(&a, -dbl) != MP_OKAY) {
         printf("\nmp_set_double() failed");
         goto LBL_ERR;
      }
      if (-dbl != mp_get_double(&a)) {
         printf("\nmp_get_double() bad result!");
         goto LBL_ERR;
      }
   }
#endif

   mp_clear_multi(&a, &b, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, NULL);
   return EXIT_FAILURE;

}

static int test_mp_get_int(void)
{
   unsigned long t;
   int i;

   mp_int a, b;
   if (mp_init_multi(&a, &b, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   for (i = 0; i < 1000; ++i) {
      t = (unsigned long)(rand() * rand() + 1) & 0xFFFFFFFFuL;
      mp_set_int(&a, t);
      if (t != mp_get_int(&a)) {
         printf("\nmp_get_int() bad result!");
         goto LBL_ERR;
      }
   }
   mp_set_int(&a, 0uL);
   if (mp_get_int(&a) != 0) {
      printf("\nmp_get_int() bad result!");
      goto LBL_ERR;
   }
   mp_set_int(&a, 0xFFFFFFFFuL);
   if (mp_get_int(&a) != 0xFFFFFFFFuL) {
      printf("\nmp_get_int() bad result!");
      goto LBL_ERR;
   }

   mp_clear_multi(&a, &b, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, NULL);
   return EXIT_FAILURE;
}

static int test_mp_get_long(void)
{
   unsigned long s, t;
   int i;

   mp_int a, b;
   if (mp_init_multi(&a, &b, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   for (i = 0; i < ((int)(sizeof(unsigned long)*CHAR_BIT) - 1); ++i) {
      t = (1ULL << (i+1)) - 1;
      if (!t)
         t = -1;
      printf(" t = 0x%lx i = %d\r", t, i);
      do {
         if (mp_set_long(&a, t) != MP_OKAY) {
            printf("\nmp_set_long() error!");
            goto LBL_ERR;
         }
         s = mp_get_long(&a);
         if (s != t) {
            printf("\nmp_get_long() bad result! 0x%lx != 0x%lx", s, t);
            goto LBL_ERR;
         }
         t <<= 1;
      } while (t != 0uL);
   }

   mp_clear_multi(&a, &b, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, NULL);
   return EXIT_FAILURE;
}

static int test_mp_get_long_long(void)
{
   unsigned long long q, r;
   int i;

   mp_int a, b;
   if (mp_init_multi(&a, &b, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   for (i = 0; i < ((int)(sizeof(unsigned long long)*CHAR_BIT) - 1); ++i) {
      r = (1ULL << (i+1)) - 1;
      if (!r)
         r = -1;
      printf(" r = 0x%llx i = %d\r", r, i);
      do {
         if (mp_set_long_long(&a, r) != MP_OKAY) {
            printf("\nmp_set_long_long() error!");
            goto LBL_ERR;
         }
         q = mp_get_long_long(&a);
         if (q != r) {
            printf("\nmp_get_long_long() bad result! 0x%llx != 0x%llx", q, r);
            goto LBL_ERR;
         }
         r <<= 1;
      } while (r != 0uLL);
   }

   mp_clear_multi(&a, &b, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, NULL);
   return EXIT_FAILURE;

}

static int test_mp_sqrt(void)
{
   int i, n;

   mp_int a, b, c, d;
   if (mp_init_multi(&a, &b, &c, &d, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   for (i = 0; i < 1000; ++i) {
      printf("%6d\r", i);
      fflush(stdout);
      n = (rand() & 15) + 1;
      mp_rand(&a, n);
      if (mp_sqrt(&a, &b) != MP_OKAY) {
         printf("\nmp_sqrt() error!");
         goto LBL_ERR;
      }
      mp_n_root_ex(&a, 2uL, &c, 0);
      mp_n_root_ex(&a, 2uL, &d, 1);
      if (mp_cmp_mag(&c, &d) != MP_EQ) {
         printf("\nmp_n_root_ex() bad result!");
         goto LBL_ERR;
      }
      if (mp_cmp_mag(&b, &c) != MP_EQ) {
         printf("mp_sqrt() bad result!\n");
         goto LBL_ERR;
      }
   }

   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_FAILURE;
}

static int test_mp_is_square(void)
{
   int i, n;

   mp_int a, b;
   if (mp_init_multi(&a, &b, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   for (i = 0; i < 1000; ++i) {
      printf("%6d\r", i);
      fflush(stdout);

      /* test mp_is_square false negatives */
      n = (rand() & 7) + 1;
      mp_rand(&a, n);
      mp_sqr(&a, &a);
      if (mp_is_square(&a, &n) != MP_OKAY) {
         printf("\nfn:mp_is_square() error!");
         goto LBL_ERR;
      }
      if (n == 0) {
         printf("\nfn:mp_is_square() bad result!");
         goto LBL_ERR;
      }

      /* test for false positives */
      mp_add_d(&a, 1uL, &a);
      if (mp_is_square(&a, &n) != MP_OKAY) {
         printf("\nfp:mp_is_square() error!");
         goto LBL_ERR;
      }
      if (n == 1) {
         printf("\nfp:mp_is_square() bad result!");
         goto LBL_ERR;
      }

   }
   printf("\n\n");

   mp_clear_multi(&a, &b, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, NULL);
   return EXIT_FAILURE;
}

static int test_mp_sqrtmod_prime(void)
{
   struct mp_sqrtmod_prime_st {
      unsigned long p;
      unsigned long n;
      mp_digit r;
   };

   static struct mp_sqrtmod_prime_st sqrtmod_prime[] = {
      { 5, 14, 3 },
      { 7, 9, 4 },
      { 113, 2, 62 }
   };
   int i;

   mp_int a, b, c;
   if (mp_init_multi(&a, &b, &c, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   /* r^2 = n (mod p) */
   for (i = 0; i < (int)(sizeof(sqrtmod_prime)/sizeof(sqrtmod_prime[0])); ++i) {
      mp_set_int(&a, sqrtmod_prime[i].p);
      mp_set_int(&b, sqrtmod_prime[i].n);
      if (mp_sqrtmod_prime(&b, &a, &c) != MP_OKAY) {
         printf("Failed executing %d. mp_sqrtmod_prime\n", (i+1));
         goto LBL_ERR;
      }
      if (mp_cmp_d(&c, sqrtmod_prime[i].r) != MP_EQ) {
         printf("Failed %d. trivial mp_sqrtmod_prime\n", (i+1));
         ndraw(&c, "r");
         goto LBL_ERR;
      }
   }

   mp_clear_multi(&a, &b, &c, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, NULL);
   return EXIT_FAILURE;
}

#if defined(LTM_DEMO_REAL_RAND) && !defined(_WIN32)
static FILE *fd_urandom = 0;
#endif

static int myrng(unsigned char *dst, int len, void *dat)
{
   int x;
   (void)dat;
#if defined(LTM_DEMO_REAL_RAND) && !defined(_WIN32)
   if (!fd_urandom) {
      fprintf(stderr, "\nno /dev/urandom\n");
   } else {
      return fread(dst, 1uL, len, fd_urandom);
   }
#endif
   for (x = 0; x < len;) {
      unsigned int r = (unsigned int)rand();
      do {
         dst[x++] = r & 0xFFu;
         r >>= 8;
      } while ((r != 0u) && (x < len));
   }
   return len;
}

static int test_mp_prime_random_ex(void)
{
   int ix, err;

   mp_int a, b;
   if (mp_init_multi(&a, &b, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   /* test for size */
   for (ix = 10; ix < 128; ix++) {
      printf("Testing (not safe-prime): %9d bits    \r", ix);
      fflush(stdout);
      err = mp_prime_random_ex(&a, 8, ix,
                               (rand() & 1) ? 0 : LTM_PRIME_2MSB_ON, myrng,
                               NULL);
      if (err != MP_OKAY) {
         printf("\nfailed with error: %s\n", mp_error_to_string(err));
         goto LBL_ERR;
      }
      if (mp_count_bits(&a) != ix) {
         printf("Prime is %d not %d bits!!!\n", mp_count_bits(&a), ix);
         goto LBL_ERR;
      }
   }
   printf("\n");

   mp_clear_multi(&a, &b, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, NULL);
   return EXIT_FAILURE;
}

static int test_mp_prime_is_prime(void)
{
   int ix, err, cnt;

   mp_int a, b;
   if (mp_init_multi(&a, &b, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   /* strong Miller-Rabin pseudoprime to the first 200 primes (F. Arnault) */
   puts("Testing mp_prime_is_prime() with Arnault's pseudoprime  803...901 \n");
   mp_read_radix(&a,
                 "91xLNF3roobhzgTzoFIG6P13ZqhOVYSN60Fa7Cj2jVR1g0k89zdahO9/kAiRprpfO1VAp1aBHucLFV/qLKLFb+zonV7R2Vxp1K13ClwUXStpV0oxTNQVjwybmFb5NBEHImZ6V7P6+udRJuH8VbMEnS0H8/pSqQrg82OoQQ2fPpAk6G1hkjqoCv5s/Yr",
                 64);
   mp_prime_is_prime(&a, 8, &cnt);
   if (cnt == MP_YES) {
      printf("Arnault's pseudoprime is not prime but mp_prime_is_prime says it is.\n");
      goto LBL_ERR;
   }
   /* About the same size as Arnault's pseudoprime */
   puts("Testing mp_prime_is_prime() with certified prime 2^1119 + 53\n");
   mp_set(&a, 1uL);
   mp_mul_2d(&a,1119,&a);
   mp_add_d(&a, 53uL, &a);
   err = mp_prime_is_prime(&a, 8, &cnt);
   /* small problem */
   if (err != MP_OKAY) {
      printf("\nfailed with error: %s\n", mp_error_to_string(err));
   }
   /* large problem */
   if (cnt == MP_NO) {
      printf("A certified prime is a prime but mp_prime_is_prime says it is not.\n");
   }
   if ((err != MP_OKAY) || (cnt == MP_NO)) {
      printf("prime tested was: ");
      mp_fwrite(&a,16,stdout);
      putchar('\n');
      goto LBL_ERR;
   }
   for (ix = 16; ix < 128; ix++) {
      printf("Testing (    safe-prime): %9d bits    \r", ix);
      fflush(stdout);
      err = mp_prime_random_ex(
               &a, 8, ix, ((rand() & 1) ? 0 : LTM_PRIME_2MSB_ON) | LTM_PRIME_SAFE,
               myrng, NULL);
      if (err != MP_OKAY) {
         printf("\nfailed with error: %s\n", mp_error_to_string(err));
         goto LBL_ERR;
      }
      if (mp_count_bits(&a) != ix) {
         printf("Prime is %d not %d bits!!!\n", mp_count_bits(&a), ix);
         goto LBL_ERR;
      }
      /* let's see if it's really a safe prime */
      mp_sub_d(&a, 1uL, &b);
      mp_div_2(&b, &b);
      err = mp_prime_is_prime(&b, 8, &cnt);
      /* small problem */
      if (err != MP_OKAY) {
         printf("\nfailed with error: %s\n", mp_error_to_string(err));
      }
      /* large problem */
      if (cnt == MP_NO) {
         printf("\nsub is not prime!\n");
      }
      if ((err != MP_OKAY) || (cnt == MP_NO)) {
         printf("prime tested was: ");
         mp_fwrite(&a,16,stdout);
         putchar('\n');
         printf("sub tested was: ");
         mp_fwrite(&b,16,stdout);
         putchar('\n');
         goto LBL_ERR;
      }

   }
   /* Check regarding problem #143 */
#ifndef MP_8BIT
   mp_read_radix(&a,
                 "FFFFFFFFFFFFFFFFC90FDAA22168C234C4C6628B80DC1CD129024E088A67CC74020BBEA63B139B22514A08798E3404DDEF9519B3CD3A431B302B0A6DF25F14374FE1356D6D51C245E485B576625E7EC6F44C42E9A63A3620FFFFFFFFFFFFFFFF",
                 16);
   err = mp_prime_strong_lucas_selfridge(&a, &cnt);
   /* small problem */
   if (err != MP_OKAY) {
      printf("\nmp_prime_strong_lucas_selfridge failed with error: %s\n", mp_error_to_string(err));
   }
   /* large problem */
   if (cnt == MP_NO) {
      printf("\n\nissue #143 - mp_prime_strong_lucas_selfridge FAILED!\n");
   }
   if ((err != MP_OKAY) || (cnt == MP_NO)) {
      printf("prime tested was: ");
      mp_fwrite(&a,16,stdout);
      putchar('\n');
      goto LBL_ERR;
   }
#endif

   printf("\n\n");

   mp_clear_multi(&a, &b, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, NULL);
   return EXIT_FAILURE;

}

static int test_mp_montgomery_reduce(void)
{
   mp_digit mp;
   int ix, i, n;
   char buf[4096];

   mp_int a, b, c, d, e;
   if (mp_init_multi(&a, &b, &c, &d, &e, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   /* test montgomery */
   for (i = 1; i <= 10; i++) {
      if (i == 10)
         i = 1000;
      printf(" digit size: %2d\r", i);
      fflush(stdout);
      for (n = 0; n < 1000; n++) {
         mp_rand(&a, i);
         a.dp[0] |= 1;

         /* let's see if R is right */
         mp_montgomery_calc_normalization(&b, &a);
         mp_montgomery_setup(&a, &mp);

         /* now test a random reduction */
         for (ix = 0; ix < 100; ix++) {
            mp_rand(&c, 1 + abs(rand()) % (2*i));
            mp_copy(&c, &d);
            mp_copy(&c, &e);

            mp_mod(&d, &a, &d);
            mp_montgomery_reduce(&c, &a, mp);
            mp_mulmod(&c, &b, &a, &c);

            if (mp_cmp(&c, &d) != MP_EQ) {
/* *INDENT-OFF* */
               printf("d = e mod a, c = e MOD a\n");
               mp_todecimal(&a, buf); printf("a = %s\n", buf);
               mp_todecimal(&e, buf); printf("e = %s\n", buf);
               mp_todecimal(&d, buf); printf("d = %s\n", buf);
               mp_todecimal(&c, buf); printf("c = %s\n", buf);
               printf("compare no compare!\n"); goto LBL_ERR;
/* *INDENT-ON* */
            }
            /* only one big montgomery reduction */
            if (i > 10) {
               n = 1000;
               ix = 100;
            }
         }
      }
   }

   printf("\n\n");

   mp_clear_multi(&a, &b, &c, &d, &e, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, &d, &e, NULL);
   return EXIT_FAILURE;

}

static int test_mp_read_radix(void)
{
   char buf[4096];

   mp_int a;
   if (mp_init_multi(&a, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   mp_read_radix(&a, "123456", 10);
   mp_toradix_n(&a, buf, 10, 3);
   printf("a == %s\n", buf);
   mp_toradix_n(&a, buf, 10, 4);
   printf("a == %s\n", buf);
   mp_toradix_n(&a, buf, 10, 30);
   printf("a == %s\n", buf);

#if 0
   for (;;) {
      fgets(buf, sizeof(buf), stdin);
      mp_read_radix(&a, buf, 10);
      mp_prime_next_prime(&a, 5, 1);
      mp_toradix(&a, buf, 10);
      printf("%s, %lu\n", buf, a.dp[0] & 3);
   }
#endif

   mp_clear_multi(&a, NULL);
   return EXIT_SUCCESS;
}

static int test_mp_cnt_lsb(void)
{
   int ix;

   mp_int a, b;
   if (mp_init_multi(&a, &b, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   mp_set(&a, 1uL);
   for (ix = 0; ix < 1024; ix++) {
      if (mp_cnt_lsb(&a) != ix) {
         printf("Failed at %d, %d\n", ix, mp_cnt_lsb(&a));
         goto LBL_ERR;
      }
      mp_mul_2(&a, &a);
   }

   mp_clear_multi(&a, &b, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, NULL);
   return EXIT_FAILURE;

}

static int test_mp_reduce_2k(void)
{
   int ix, cnt;

   mp_int a, b, c, d;
   if (mp_init_multi(&a, &b, &c, &d, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   /* test mp_reduce_2k */
   for (cnt = 3; cnt <= 128; ++cnt) {
      mp_digit tmp;

      mp_2expt(&a, cnt);
      mp_sub_d(&a, 2uL, &a);  /* a = 2**cnt - 2 */

      printf("\r %4d bits", cnt);
      printf("(%d)", mp_reduce_is_2k(&a));
      mp_reduce_2k_setup(&a, &tmp);
      printf("(%lu)", (unsigned long) tmp);
      for (ix = 0; ix < 1000; ix++) {
         if (!(ix & 127)) {
            printf(".");
            fflush(stdout);
         }
         mp_rand(&b, (cnt / DIGIT_BIT + 1) * 2);
         mp_copy(&c, &b);
         mp_mod(&c, &a, &c);
         mp_reduce_2k(&b, &a, 2uL);
         if (mp_cmp(&c, &b) != MP_EQ) {
            printf("FAILED\n");
            goto LBL_ERR;
         }
      }
   }

   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, &d, NULL);
   return EXIT_FAILURE;
}

static int test_mp_div_3(void)
{
   int cnt;

   mp_int a, b, c, d, e;
   if (mp_init_multi(&a, &b, &c, &d, &e, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   /* test mp_div_3  */
   mp_set(&d, 3uL);
   for (cnt = 0; cnt < 10000;) {
      mp_digit r2;

      if (!(++cnt & 127)) {
         printf("%9d\r", cnt);
         fflush(stdout);
      }
      mp_rand(&a, abs(rand()) % 128 + 1);
      mp_div(&a, &d, &b, &e);
      mp_div_3(&a, &c, &r2);

      if (mp_cmp(&b, &c) || mp_cmp_d(&e, r2)) {
         printf("\nmp_div_3 => Failure\n");
         goto LBL_ERR;
      }
   }
   printf("\nPassed div_3 testing");

   mp_clear_multi(&a, &b, &c, &d, &e, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, &d, &e, NULL);
   return EXIT_FAILURE;
}

static int test_mp_dr_reduce(void)
{
   mp_digit mp;
   int cnt;
   unsigned rr;
   int ix;

   mp_int a, b, c;
   if (mp_init_multi(&a, &b, &c, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }

   /* test the DR reduction */
   for (cnt = 2; cnt < 32; cnt++) {
      printf("\r%d digit modulus", cnt);
      mp_grow(&a, cnt);
      mp_zero(&a);
      for (ix = 1; ix < cnt; ix++) {
         a.dp[ix] = MP_MASK;
      }
      a.used = cnt;
      a.dp[0] = 3;

      mp_rand(&b, cnt - 1);
      mp_copy(&b, &c);

      rr = 0;
      do {
         if (!(rr & 127)) {
            printf(".");
            fflush(stdout);
         }
         mp_sqr(&b, &b);
         mp_add_d(&b, 1uL, &b);
         mp_copy(&b, &c);

         mp_mod(&b, &a, &b);
         mp_dr_setup(&a, &mp);
         mp_dr_reduce(&c, &a, mp);

         if (mp_cmp(&b, &c) != MP_EQ) {
            printf("Failed on trial %u\n", rr);
            goto LBL_ERR;
         }
      } while (++rr < 500);
      printf(" passed");
      fflush(stdout);
   }

   mp_clear_multi(&a, &b, &c, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, &c, NULL);
   return EXIT_FAILURE;
}

static int test_mp_reduce_2k_l(void)
{
#   if LTM_DEMO_TEST_REDUCE_2K_L
   mp_int a, b;
   if (mp_init_multi(&a, &b, NULL)!= MP_OKAY) {
      return EXIT_FAILURE;
   }
   /* test the mp_reduce_2k_l code */
#      if LTM_DEMO_TEST_REDUCE_2K_L == 1
   /* first load P with 2^1024 - 0x2A434 B9FDEC95 D8F9D550 FFFFFFFF FFFFFFFF */
   mp_2expt(&a, 1024);
   mp_read_radix(&b, "2A434B9FDEC95D8F9D550FFFFFFFFFFFFFFFF", 16);
   mp_sub(&a, &b, &a);
#      elif LTM_DEMO_TEST_REDUCE_2K_L == 2
   /*  p = 2^2048 - 0x1 00000000 00000000 00000000 00000000 4945DDBF 8EA2A91D 5776399B B83E188F  */
   mp_2expt(&a, 2048);
   mp_read_radix(&b,
                 "1000000000000000000000000000000004945DDBF8EA2A91D5776399BB83E188F",
                 16);
   mp_sub(&a, &b, &a);
#      else
#         error oops
#      endif

   mp_todecimal(&a, buf);
   printf("\n\np==%s\n", buf);
   /* now mp_reduce_is_2k_l() should return */
   if (mp_reduce_is_2k_l(&a) != 1) {
      printf("mp_reduce_is_2k_l() return 0, should be 1\n");
      goto LBL_ERR;
   }
   mp_reduce_2k_setup_l(&a, &d);
   /* now do a million square+1 to see if it varies */
   mp_rand(&b, 64);
   mp_mod(&b, &a, &b);
   mp_copy(&b, &c);
   printf("Testing: mp_reduce_2k_l...");
   fflush(stdout);
   for (cnt = 0; cnt < (int)(1UL << 20); cnt++) {
      mp_sqr(&b, &b);
      mp_add_d(&b, 1uL, &b);
      mp_reduce_2k_l(&b, &a, &d);
      mp_sqr(&c, &c);
      mp_add_d(&c, 1uL, &c);
      mp_mod(&c, &a, &c);
      if (mp_cmp(&b, &c) != MP_EQ) {
         printf("mp_reduce_2k_l() failed at step %d\n", cnt);
         mp_tohex(&b, buf);
         printf("b == %s\n", buf);
         mp_tohex(&c, buf);
         printf("c == %s\n", buf);
         goto LBL_ERR;
      }
   }

   mp_clear_multi(&a, &b, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &b, NULL);
   return EXIT_FAILURE;
#else
   return EXIT_SUCCESS;
#   endif /* LTM_DEMO_TEST_REDUCE_2K_L */
}

#ifdef LTM_USE_EXTRA_FUNCTIONS
#if ( (defined LTM_USE_FASTER_VERSIONS) || (defined LTM_USE_FASTER_RADIX_SIZE) )
/* stripped down version of mp_radix_size. The faster version can be off by up to +3  */
static int s_rs(const mp_int *a, int radix, int *size)
{
   int     res, digs = 0;
   mp_int  t;
   mp_digit d;
   *size = 0;
   if (mp_iszero(a) == MP_YES) {
      *size = 2;
      return MP_OKAY;
   }
   if (radix == 2) {
      *size = mp_count_bits(a) + 1;
      return MP_OKAY;
   }
   if ((res = mp_init_copy(&t, a)) != MP_OKAY) {
      return res;
   }
   t.sign = MP_ZPOS;
   while (mp_iszero(&t) == MP_NO) {
      if ((res = mp_div_d(&t, (mp_digit)radix, &t, &d)) != MP_OKAY) {
         mp_clear(&t);
         return res;
      }
      ++digs;
   }
   mp_clear(&t);
   *size = digs + 1;
   return MP_OKAY;
}
#endif

static int test_mp_ilogb(void)
{
   mp_int a, lb;
   mp_digit d, base;
   int size;

   mp_init_multi(&a, &lb, NULL);

   /*
     base   a    result
      0     x    MP_VAL
      1     x    MP_VAL
   */
   mp_set(&a, 42uL);
   base = 0uL;
   if (mp_ilogb(&a, base, &lb) != MP_VAL) {
      goto LBL_ERR;
   }
   base = 1uL;
   if (mp_ilogb(&a, base, &lb) != MP_VAL) {
      goto LBL_ERR;
   }
   /*
     base   a    result
      2     0    MP_VAL
      2     1    0
      2     2    1
      2     3    1
   */
   base = 2uL;
   mp_zero(&a);
   if (mp_ilogb(&a, base, &lb) != MP_VAL) {
      goto LBL_ERR;
   }

   for (d = 1; d < 4; d++) {
      mp_set(&a, d);
      if (mp_ilogb(&a, base, &lb) != MP_OKAY) {
         goto LBL_ERR;
      }
      if (mp_cmp_d(&lb, (d == 1)?0uL:1uL) != MP_EQ) {
         goto LBL_ERR;
      }
   }
   /*
    base   a    result
     3     0    MP_VAL
     3     1    0
     3     2    0
     3     3    1
   */
   base = 3uL;
   mp_zero(&a);
   if (mp_ilogb(&a, base, &lb) != MP_VAL) {
      goto LBL_ERR;
   }
   for (d = 1; d < 4; d++) {
      mp_set(&a, d);
      if (mp_ilogb(&a, base, &lb) != MP_OKAY) {
         goto LBL_ERR;
      }
      if (mp_cmp_d(&lb, (d < base)?0uL:1uL) != MP_EQ) {
         goto LBL_ERR;
      }
   }

   /*
     bases 2..64 with "a" a random large constant.
     The range of bases tested allows to check with
     radix_size.
   */
   mp_rand(&a, 10);
   for (base = 2uL; base < 65uL; base++) {
      if (mp_ilogb(&a, base, &lb) != MP_OKAY) {
         goto LBL_ERR;
      }
#if ( (defined LTM_USE_FASTER_VERSIONS) || (defined LTM_USE_FASTER_RADIX_SIZE) )
      if (s_rs(&a,(int)base, &size) != MP_OKAY) {
         goto LBL_ERR;
      }
#else
      if (mp_radix_size(&a,(int)base, &size) != MP_OKAY) {
         goto LBL_ERR;
      }
#endif
      /* radix_size includes the memory needed for '\0', too*/
      size -= 2;
      if (mp_cmp_d(&lb, size) != MP_EQ) {
         goto LBL_ERR;
      }
   }

   /*Test upper edgecase with base MP_MASK and number (MP_MASK/2)*MP_MASK^10  */
   mp_set(&a, MP_MASK);
   if (mp_expt_d(&a, 10uL, &a) != MP_OKAY) {
      goto LBL_ERR;
   }
   if (mp_add_d(&a, (MP_MASK>>1), &a) != MP_OKAY) {
      goto LBL_ERR;
   }
   if (mp_ilogb(&a, MP_MASK, &lb) != MP_OKAY) {
      goto LBL_ERR;
   }
   if (mp_cmp_d(&lb, 10uL) != MP_EQ) {
      goto LBL_ERR;
   }

   mp_clear_multi(&a, &lb, NULL);
   return EXIT_SUCCESS;
LBL_ERR:
   mp_clear_multi(&a, &lb, NULL);
   return EXIT_FAILURE;
}
static int test_mp_is_small_prime(void)
{
   mp_sieve *base = NULL;
   mp_sieve *segment = NULL;
   LTM_SIEVE_UINT single_segment_a = 0;
   int e;
   int i, test_size;

   LTM_SIEVE_UINT to_test[] = {
      52, 137, 153, 179, 6, 153, 53, 132, 150, 65,
      27414, 36339, 36155, 11067, 52060, 5741,
      29755, 2698, 52572, 13053, 9375, 47241
#ifndef MP_8BIT
      ,39626, 207423, 128857, 37419, 141696, 189465,
      41503, 127370, 91673, 8473, 479142414, 465566339,
      961126169, 1057886067, 1222702060, 1017450741,
      1019879755, 72282698, 2048787577, 2058368053
#endif
   };
   LTM_SIEVE_UINT tested[] = {
      0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0
#ifndef MP_8BIT
      ,0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0
#endif
   };
   LTM_SIEVE_UINT result;

   if ((e = mp_sieve_init(base)) != MP_OKAY) {
      fprintf(stderr,"mp_sieve_init(base) failed: \"%s\"\n",mp_error_to_string(e));
      exit(EXIT_FAILURE);
   }
   if ((e = mp_sieve_init(segment)) != MP_OKAY) {
      fprintf(stderr,"mp_sieve_init(segment) failed: \"%s\"\n",mp_error_to_string(e));
      goto LTM_ERR;
   }

   test_size = (int)(sizeof(to_test)/sizeof(LTM_SIEVE_UINT));

   for (i = 0; i < test_size; i++) {
      if ((e = mp_is_small_prime(to_test[i], &result, &base, &segment, &single_segment_a)) != MP_OKAY) {
         fprintf(stderr,"mp_is_small_prime failed: \"%s\"\n",mp_error_to_string(e));
         goto LTM_ERR;
      }
      if (result != tested[i]) {
         fprintf(stderr,"mp_is_small_prime failed for %u. Said %u but is %u \n",
                 (unsigned int)to_test[i], (unsigned int)result, (unsigned int)tested[i]);
         goto LTM_ERR;
      }
   }

   mp_sieve_clear(base);
   mp_sieve_clear(segment);
   return EXIT_SUCCESS;
LTM_ERR:
   mp_sieve_clear(base);
   mp_sieve_clear(segment);
   return EXIT_FAILURE;
}
static int test_mp_next_small_prime(void)
{
   mp_sieve *base = NULL;
   mp_sieve *segment = NULL;
   LTM_SIEVE_UINT single_segment_a = 0;
   LTM_SIEVE_UINT ret;
   int e;
   int i, test_size;

   LTM_SIEVE_UINT to_test[] = {
      52, 137, 153, 179, 6, 153, 53, 132, 150, 65,
      27414, 36339, 36155, 11067, 52060, 5741,
      29755, 2698, 52572, 13053, 9375, 47241
#ifndef MP_8BIT
      ,39626, 207423, 128857, 37419, 141696, 189465,
      41503, 127370, 91673, 8473, 479142414, 465566339,
      961126169, 1057886067, 1222702060, 1017450741,
      1019879755, 72282698, 2048787577, 2058368053
#endif
   };
   LTM_SIEVE_UINT tested[] = {
      53, 137, 157, 179, 7, 157, 53, 137, 151, 67,
      27427, 36341, 36161, 11069, 52067, 5741,
      29759, 2699, 52579, 13063, 9377, 47251
#ifndef MP_8BIT
      ,39631, 207433, 128857, 37423, 141697, 189467,
      41507, 127373, 91673, 8501, 479142427, 465566393,
      961126169, 1057886083, 1222702081, 1017450823,
      1019879761, 72282701, 2048787577, 2058368113
#endif
   };

   if ((e = mp_sieve_init(base)) != MP_OKAY) {
      fprintf(stderr,"mp_sieve_init(base) failed: \"%s\"\n",mp_error_to_string(e));
      return EXIT_FAILURE;
   }
   if ((e = mp_sieve_init(segment)) != MP_OKAY) {
      fprintf(stderr,"mp_sieve_init(segment) failed: \"%s\"\n",mp_error_to_string(e));
      goto LTM_ERR_1;
   }

   test_size = (int)(sizeof(to_test)/sizeof(LTM_SIEVE_UINT));

   for (i = 0; i < test_size; i++) {
      if ((e = mp_next_small_prime(to_test[i], &ret, &base,
                                   &segment, &single_segment_a)) != MP_OKAY) {
         fprintf(stderr,"mp_next_small_prime failed with \"%s\" at index %d\n",
                 mp_error_to_string(e), i);
         goto LTM_ERR;
      }
      if (ret != tested[i]) {
         fprintf(stderr,"mp_next_small_prime failed for %u. Said %u but is %u \n",
                 (unsigned int)to_test[i], (unsigned int)ret, (unsigned int)tested[i]);
         goto LTM_ERR;
      }
   }

   mp_sieve_clear(base);
   mp_sieve_clear(segment);
   return EXIT_SUCCESS;
LTM_ERR:
   mp_sieve_clear(segment);
LTM_ERR_1:
   mp_sieve_clear(base);
   return EXIT_FAILURE;
}

static int test_mp_prec_small_prime(void)
{
   mp_sieve *base = NULL;
   mp_sieve *segment = NULL;
   LTM_SIEVE_UINT single_segment_a = 0;
   LTM_SIEVE_UINT ret;
   int e;
   int i, test_size;

   LTM_SIEVE_UINT to_test[] = {
      52, 137, 153, 179, 6, 153, 53, 132, 150, 65,
      27414, 36339, 36155, 11067, 52060, 5741,
      29755, 2698, 52572, 13053, 9375, 47241
#ifndef MP_8BIT
      ,39626, 207423, 128857, 37419, 141696, 189465,
      41503, 127370, 91673, 8473, 479142414, 465566339,
      961126169, 1057886067, 1222702060, 1017450741,
      1019879755, 72282698, 2048787577, 2058368053
#endif
   };
   LTM_SIEVE_UINT tested[] = {
      47, 137, 151, 179, 5, 151, 53, 131, 149, 61,
      27409, 36319, 36151, 11059, 52057, 5741,
      29753, 2693, 52571, 13049, 9371, 47237
#ifndef MP_8BIT
      ,39623, 207409, 128857, 37409, 141689, 189463,
      41491, 127363, 91673, 8467, 479142413, 465566323,
      961126169, 1057886029, 1222702051, 1017450739,
      1019879717, 72282697, 2048787577, 2058368051
#endif
   };

   if ((e = mp_sieve_init(base)) != MP_OKAY) {
      fprintf(stderr,"mp_sieve_init(base) failed: \"%s\"\n",mp_error_to_string(e));
      return EXIT_FAILURE;
   }
   if ((e = mp_sieve_init(segment)) != MP_OKAY) {
      fprintf(stderr,"mp_sieve_init(segment) failed: \"%s\"\n",mp_error_to_string(e));
      goto LTM_ERR_1;
   }

   test_size = (int)(sizeof(to_test)/sizeof(LTM_SIEVE_UINT));

   for (i = 0; i < test_size; i++) {
      if ((e = mp_prec_small_prime(to_test[i], &ret, &base,
                                   &segment, &single_segment_a)) != MP_OKAY) {
         fprintf(stderr,"mp_prec_small_prime failed with \"%s\" at index %d\n",
                 mp_error_to_string(e), i);
         goto LTM_ERR;
      }
      if (ret != tested[i]) {
         fprintf(stderr,"mp_prec_small_prime failed for %u. Said %u but is %u \n",
                 (unsigned int)to_test[i], (unsigned int)ret, (unsigned int)tested[i]);
         goto LTM_ERR;
      }
   }

   mp_sieve_clear(base);
   mp_sieve_clear(segment);
   return EXIT_SUCCESS;
LTM_ERR:
   mp_sieve_clear(segment);
LTM_ERR_1:
   mp_sieve_clear(base);
   return EXIT_FAILURE;
}

static int test_mp_small_prime_array(void)
{
   mp_sieve *base = NULL;
   mp_sieve *segment = NULL;
   LTM_SIEVE_UINT single_segment_a = 0;
   LTM_SIEVE_UINT *prime_array, array_size, j;
   int e;
   int i, test_size;
   mp_int total, t;

   /* For simplicity: just check the primesums */
   LTM_SIEVE_UINT to_test_start[] = {
      124, 147, 185, 291, 311, 329
#ifndef MP_8BIT
      , 7901, 19489, 20524, 47371, 50741, 59233
#endif
   };
   LTM_SIEVE_UINT to_test_end[] = {
      479, 488, 661, 671, 703, 833
#ifndef MP_8BIT
      , 66339, 68053, 72241, 81055, 82698, 86067
#endif
   };
   LTM_SIEVE_UINT tested[] = {
      18466, 18419, 33441, 28906, 31731, 45475
#ifndef MP_8BIT
      , 203371931, 197821202, 221915352, 195494605, 192053247, 172585395
#endif
   };
   /* TODO: add a pair or two that overlap segments */

   if ((e = mp_sieve_init(base)) != MP_OKAY) {
      fprintf(stderr, "mp_sieve_init(base) failed: \"%s\"\n", mp_error_to_string(e));
      return EXIT_FAILURE;
   }
   if ((e = mp_sieve_init(segment)) != MP_OKAY) {
      fprintf(stderr, "mp_sieve_init(segment) failed: \"%s\"\n", mp_error_to_string(e));
      goto LTM_ERR_2;
   }

   if ((e = mp_init_multi(&total, &t, NULL)) != MP_OKAY) {
      fprintf(stderr, "mp_init_multi(segment): \"%s\"\n", mp_error_to_string(e));
      goto LTM_ERR_1;
   }

   test_size = (int)(sizeof(to_test_start) / sizeof(LTM_SIEVE_UINT));

   for (i = 0; i < test_size; i++) {
      mp_zero(&total);
      if ((e = mp_small_prime_array(to_test_start[i], to_test_end[i], &prime_array,
                                    &array_size, &base,
                                    &segment, &single_segment_a)) != MP_OKAY) {
         fprintf(stderr, "mp_small_prime_array failed with \"%s\" for %u %u\n",
                 mp_error_to_string(e), to_test_start[i], to_test_end[i]);
         goto LTM_ERR;
      }
      for (j = 0; j < array_size; j++) {
#ifdef MP_64BIT
         if ((e = mp_add_d(&total, (mp_digit)prime_array[j], &total)) != MP_OKAY) {
            fprintf(stderr,"mp_add_d failed: \"%s\"\n",mp_error_to_string(e));
            goto LTM_ERR;
         }
#else
         if ((e = mp_set_long(&t,(mp_word)prime_array[j])) != MP_OKAY) {
            fprintf(stderr,"mp_set_long (1) failed: \"%s\"\n",mp_error_to_string(e));
            goto LTM_ERR;
         }
         if ((e = mp_add(&total, &t, &total)) != MP_OKAY) {
            fprintf(stderr,"mp_add failed: \"%s\"\n",mp_error_to_string(e));
            goto LTM_ERR;
         }
#endif
      }
#ifdef MP_64BIT
      if (mp_cmp_d(&total, (mp_digit)tested[i]) != MP_EQ) {
         fprintf(stderr,
                 "mp_small_prime_array primesum failed for %u %u == %u not ",
                 to_test_start[i], to_test_end[i], tested[i]);
         goto LTM_ERR;
      }
#else
      if ((e = mp_set_long(&t,(mp_word)tested[i])) != MP_OKAY) {
          fprintf(stderr,"mp_set_long (1) failed: \"%s\"\n",mp_error_to_string(e));
          goto LTM_ERR;
       }
      if (mp_cmp(&total, &t) != MP_EQ) {
         fprintf(stderr,
                 "mp_small_prime_array (8-bit) primesum failed for %u %u == %u not ",
                 to_test_start[i], to_test_end[i], tested[i]);
         mp_fwrite(&total, 10,stdout);
         putchar('\n');
         goto LTM_ERR;
      }
#endif
   }

   mp_clear_multi(&total, &t, NULL);
   mp_sieve_clear(base);
   mp_sieve_clear(segment);
   return EXIT_SUCCESS;
LTM_ERR:
   mp_clear_multi(&total, &t, NULL);
LTM_ERR_1:
   mp_sieve_clear(segment);
LTM_ERR_2:
   mp_sieve_clear(base);
   return EXIT_FAILURE;
}

#endif

int unit_tests(void)
{
   static const struct {
      const char *name;
      int (*fn)(void);
   } test[] = {
#define T(n) { #n, test_##n }
            T(trivial_stuff),
            T(mp_cnt_lsb),
            T(mp_complement),
            T(mp_div_3),
            T(mp_dr_reduce),
            T(mp_get_int),
            T(mp_get_long),
            T(mp_get_long_long),
            T(mp_invmod),
            T(mp_is_square),
            T(mp_jacobi),
            T(mp_kronecker),
            T(mp_montgomery_reduce),
            T(mp_prime_is_prime),
            T(mp_prime_random_ex),
            T(mp_read_radix),
            T(mp_reduce_2k),
            T(mp_reduce_2k_l),
            T(mp_set_double),
            T(mp_sqrt),
            T(mp_sqrtmod_prime),
            T(mp_tc_and),
            T(mp_tc_div_2d),
            T(mp_tc_or),
            T(mp_tc_xor)
#ifdef LTM_USE_EXTRA_FUNCTIONS
     ,T(mp_ilogb),
      T(mp_is_small_prime),
      T(mp_next_small_prime),
      T(mp_prec_small_prime),
      T(mp_small_prime_array)
#endif
#undef T
   };
   unsigned long i;
   int res = EXIT_SUCCESS;

#if defined(LTM_DEMO_REAL_RAND) && !defined(_WIN32)
   fd_urandom = fopen("/dev/urandom", "r");
   if (!fd_urandom) {
      fprintf(stderr, "\ncould not open /dev/urandom\n");
   }
#endif

   for (i = 0; i < sizeof(test) / sizeof(test[0]); ++i) {
      printf("TEST %s ", test[i].name);
      if (test[i].fn() != EXIT_SUCCESS) {
         puts("FAIL");
         res = EXIT_FAILURE;
         break;
      }
      puts("PASS");
   }

#if defined(LTM_DEMO_REAL_RAND) && !defined(_WIN32)
   if (fd_urandom) {
      fclose(fd_urandom);
   }
#endif
   return res;
}

