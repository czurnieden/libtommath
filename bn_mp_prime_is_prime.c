#include <tommath.h>
#ifdef BN_MP_PRIME_IS_PRIME_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis
 *
 * LibTomMath is a library that provides multiple-precision
 * integer arithmetic as well as number theoretic functionality.
 *
 * The library was designed directly after the MPI library by
 * Michael Fromberger but has been written from scratch with
 * additional optimizations in place.
 *
 * The library is free for all purposes without any express
 * guarantee it works.
 *
 * Tom St Denis, tomstdenis@gmail.com, http://libtom.org
 */

/* performs a variable number of rounds of Miller-Rabin
 *
 * Probability of error after t rounds is no more than

 *
 * Sets result to 1 if probably prime, 0 otherwise
 */


/*
 * Port from the function iLucasSelfridge(mpz_t mpzN) in trn.c, which is:
 *
 * Freeware copyright (c) 2012 Thomas R. Nicely <http://www.trnicely.net>.
 * Released into the public domain by the author, who disclaims any legal
 * liability arising from its use.
 *
 * Ported to Calc by Christoph Zurnieden.
 */
#if 0
static unsigned long uLDmax = 0;
static int lucasselfridge(mp_int zN)
{
    // the "long" is assumed to have 32 bit?
    int cmpz, iP, iJ, iSign, ret;
    long lDabs, lD, lQ;
    unsigned long ulNbits, ul, ulGCD;
    mp_int zU, zV, zNplus1, zU2m, zV2m, zQm, z2Qm, zT1, zT2, zT3, zT4, zD;

    mp_init_multi(&zU, &zV, &zNplus1, &zU2m, &zV2m, &zQm, &z2Qm, &zT1, &zT2, &zT3, &zT4, &zD, NULL);

    // May I point you to our brand new product "mp_zero_multi()"? No? Well...
    mp_zero(&zU);
    mp_zero(&zV);
    mp_zero(&zNplus1);
    mp_zero(&zU2m);
    mp_zero(&zV2m);
    mp_zero(&zQm);
    mp_zero(&z2Qm);
    mp_zero(&zT1);
    mp_zero(&zT2);
    mp_zero(&zT3);
    mp_zero(&zT4);
    mp_zero(&zD);

// If you think that's an abomination you must see the original.
#undef RETURN
#define RETURN(n)           \
  {                         \
  mp_clear_multi(&zU, &zV, &zNplus1, &zU2m, &zV2m, &zQm, &z2Qm, &zT1, &zT2, &zT3, &zT4, &zD, NULL)\
  return(n);                \
  }


    cmpz = mp_cmp_d(&zN, (mp_digit) 2);
    if (cmpz == MP_LT) {
	return 0;
    }
    if (cmpz == MP_EQ) {
	return 1;
    }
    if (mp_iseven(&zN)) {
	return 0;
    }
    mp_is_square(&zN, &ret);
    if (ret != 0) {
	return 0;
    }

    lDabs = 5;
    iSign = 1;
#include <bool.h>
    while (true) {
	lD = iSign * lDabs;
	iSign = -iSign;

	ulGCD = zgcd_ui(zN, lDabs, NULL);



	if ((ulGCD > 1) && zrel_ui(zN, ulGCD) > 0) {
	    RETURN(0);
	}

	itoz(lD, &zD);

	iJ = zjacobi(zD, zN);

	if (iJ == -1) {
	    break;
	}

	lDabs += 2;
	if (lDabs > (long) uLDmax) {
	    uLDmax = lDabs;
	}
	if (lDabs > INT32_MAX - 2) {
	    fprintf(stderr,
		    "\n ERROR: D overflows signed long in Lucas-Selfridge test.");
	    fprintf(stderr, "\n N=");
	    zfprintf(stderr, zN, 0, 0);
	    fprintf(stderr, "\n |D|=%ld\n\n", lDabs);
	    exit(EXIT_FAILURE);
	}
    }

    iP = 1;
    lQ = (1 - lD) / 4;

    zadd_ui(zN, 1, &zNplus1);

    utoz(0, &zU);
    utoz(2, &zV);
    utoz(1, &zU2m);
    itoz(iP, &zV2m);
    itoz(lQ, &zQm);
    itoz(2 * lQ, &z2Qm);

    ulNbits = zhighbit(zNplus1) + 1;
    for (ul = 1; ul < ulNbits; ul++) {

	zmul(zU2m, zV2m, &zU2m);
	zmod(zU2m, zN, &zU2m, 0);
	zmul(zV2m, zV2m, &zV2m);
	zsub(zV2m, z2Qm, &zV2m);
	zmod(zV2m, zN, &zV2m, 0);

	if (zisset(zNplus1, ul)) {
	    zmul(zU2m, zV, &zT1);
	    zmul(zU, zV2m, &zT2);
	    zmul(zV2m, zV, &zT3);
	    zmul(zU2m, zU, &zT4);

	    zmuli(zT4, lD, &zT4);

	    zadd(zT1, zT2, &zU);

	    if (zisodd(zU)) {
		zadd(zU, zN, &zU);
	    }

	    /*mpz_fdiv_q_2exp(zU, zU, 1); */
	    zdivi(zU, 2, &zU);

	    zadd(zT3, zT4, &zV);

	    if (zisodd(zV)) {
		zadd(zV, zN, &zV);
	    }

	    /*mpz_fdiv_q_2exp(zV, zV, 1); */
	    zdivi(zV, 2, &zV);

	    zmod(zU, zN, &zU, 0);
	    zmod(zV, zN, &zV, 0);

	}
	if (ul < ulNbits - 1) {
	    zmul(zQm, zQm, &zQm);
	    zmod(zQm, zN, &zQm, 0);
	    zadd(zQm, zQm, &z2Qm);
	}

    }

    if (ziszero(zU)) {
	RETURN(1);
    }
    RETURN(0);
}

#endif
int mp_prime_is_prime(mp_int *a, int t, int *result)
{
   mp_int  b;
   int     ix, err, res;

   /* default to no */
   *result = MP_NO;

   /* valid value of t? */
   if (t <= 0 || t > PRIME_SIZE) {
      return MP_VAL;
   }

   /* is the input equal to one of the primes in the table? */
   for (ix = 0; ix < PRIME_SIZE; ix++) {
      if (mp_cmp_d(a, ltm_prime_tab[ix]) == MP_EQ) {
         *result = 1;
         return MP_OKAY;
      }
   }

   /* first perform trial division */
   if ((err = mp_prime_is_divisible(a, &res)) != MP_OKAY) {
      return err;
   }

   /* return if it was trivially divisible */
   if (res == MP_YES) {
      return MP_OKAY;
   }

   /* now perform the miller-rabin rounds */
   if ((err = mp_init(&b)) != MP_OKAY) {
      return err;
   }

   for (ix = 0; ix < t; ix++) {
      /* set the prime */
      mp_set(&b, ltm_prime_tab[ix]);

      if ((err = mp_prime_miller_rabin(a, &b, &res)) != MP_OKAY) {
         goto LBL_B;
      }

      if (res == MP_NO) {
         goto LBL_B;
      }
   }

   /* passed the test */
   *result = MP_YES;
LBL_B:
   mp_clear(&b);
   return err;
}
#endif

/* $Source$ */
/* $Revision$ */
/* $Date$ */
