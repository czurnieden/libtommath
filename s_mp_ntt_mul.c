#include "tommath_private.h"
#ifdef S_MP_NTT_MUL_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

/* MP_16BIT is just too small, not worth it */
#if (MP_DIGIT_BIT > 16)

/*
   TODO: adapt to LTM. ("-" nothing done, "|" partially done, "-" not done)
    | use mp_err to handle returns.
    + convert data types to mp_digit/mp_word
    | variable declarations at the beginning of the function
    + remove redundant functions
    - replace functions that are in LTM already
    - replace logical ints with bools
    + change names according to LTM rules
    + malloc -> MP_MALLOC, free -> MP_FREE
    - put individual functions into individual files?
*/
/*
  Each NTT prime p < 2^62, so:
    - residues fit in a 64-bit word
    - the product of two residues is < 2^124 and fits in a 128 bit mp_word
    - Montgomery reduction with R = 2^64 keeps every intermediate < 2^128

   This would, theoretically, fit into MP_32BIT with 28 bit limbsize, too.

   If each NTT prime is below 2^30 (assuming an unsigned 32 bit mp_digit)
    - residues fit in a 32-bit word
    - the product of two residues is < 2^60 and fits in a 64 bit mp_word

   One problem is to find good primes inside that lower limit.
   Another one: does this construct support 2^31 bit long products with three primes?
   No, it does not, it needs four.

 */

mp_digit s_mp_ntt_invmod_d(mp_digit p)
{
   mp_digit x = 1;
   int i;
   for (i = 0; i < 6; i++){
      x = x * (mp_digit)(2 - p * x);
   }
   return x;
}

/* TODO: These three (mulmod, powmod and miller-rabin) are in the branch 'all_in'  */

/* Non-Montgomery s_mp_small_mulmod/s_mp_small_powmod, used for prime test. */
mp_digit s_mp_ntt_small_mulmod(mp_digit a, mp_digit b, mp_digit m)
{
   return (mp_digit)((mp_word)a * b % m);
}

mp_digit s_mp_ntt_small_powmod(mp_digit b, mp_digit e, mp_digit m)
{
   mp_digit r = 1u;
   b %= m;
   while (e) {
      if (e & 1) {
         r = s_mp_ntt_small_mulmod(r, b, m);
      }
      b = s_mp_ntt_small_mulmod(b, b, m);
      e >>= 1;
   }
   return r;
}

/* Trial division with following Miller-Rabin with enough witnesses to check all primes < 2^64 */
bool s_mp_ntt_is_prime(mp_digit n)
{
   /* More, but all primes < 2^15 */
   /* const mp_digit bases[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}; */
   /* All primes < 2^31 */
   const mp_digit bases[7] = {2u, 325u, 9375u, 28178u, 450775u, 9780504u, 1795265022u};
   int i, j, r;
   int composite;
   mp_digit a, x, d;

   if (n < 2u) {
      return false;
   }
   if (n == 2u) {
      return true;
   }
   if ((n & 1) == 0) {
      return false;
   }
   /* Yes, my dear linter, these are magic numbers, get used to it */
   for (mp_digit p = 3u; (p < 1000u) && (p * p <= n); p+= 2u) {
      if (n % p == 0) {
         /* In case of equality n is prime. Why?
            All composites are made from smaller primes,
            so a composite p cannot be equal to n because
            smaller factors of that composite prime divide
            n and would have come earlier in the loop.
            Only if n itself is prime n == p is true.
          */
         return (n == p);
      }
   }

   d = n - 1;
   r = 0;
   while (!(d & 1)) {
      d >>= 1;
      r++;
   }
   for (i = 0; i < 7; i++) {
      a = bases[i] % n;
      if (a == 0) {
         continue;
      }
      x = s_mp_ntt_small_powmod(a, d, n);
      if (x == 1 || x == n - 1) {
         continue;
      }
      composite = 1;
      for (j = 1; j < r; j++) {
         x = s_mp_ntt_small_mulmod(x, x, n);
         if (x == n - 1) {
            composite = 0;
            break;
         }
      }
      if (composite) {
         return false;
      }
   }
   return true;
}




/*   We rely on p < 2^63, which the NTT primes (< 2^62) always satisfy:
 *   - REDC's output r is in [0, 2p) and every input a, b, s is in [0, p);
 *   - subtracting p from a value below p wraps to a value whose top bit is
 *     set, because 2^64 - p > 2^63.  The borrow is therefore just bit 63 of
 *     the wrapped difference, and "add p back if we borrowed" is one mask and
 *     one add.  The same expression would NOT be valid for a modulus close to
 *     2^64 (e.g. 2^64 - 2^32 + 1): there 2^64 - p is small and the wrapped
 *     difference does not reach bit 63.
 */

/* TODO: one does not need a profiler get that this function is the most called
         function that does anything useful. Optimizing here with soome
         assembler if the architecture supports it.
         (e.g., Intel's AVX512 madd52* (madd52hi and madd52lo) in 256-bit registers,)
 */

/* REDC: T in [0, p*R) -> T * R^{-1} mod p, in [0, p). */
mp_digit s_mp_ntt_mont_redc(const ntt_mont *m, mp_word T)
{
   mp_digit q, r, d;
   mp_word u;
   q = (mp_digit)T * m->ninv;               /* (T mod R) * (-p^{-1}) mod R */
   u = T + (mp_word)q * m->p;               /* u < 2*p*R < 2^127 */
   r = (mp_digit)(u >> NTT_BASE_SIZE);      /* r < 2*p */
   d = r - m->p;
   d = d + (m->p & ((mp_digit)0u - (d >> NTT_BASE_SIZE_M_ONE)));
   return d;
}

/* Montgomery multiply: a,b in [0,p) -> a*b*R^{-1} mod p. */
mp_digit s_mp_ntt_mont_mul(const ntt_mont *m, mp_digit a, mp_digit b)
{
   return s_mp_ntt_mont_redc(m, (mp_word)a * b);
}

/* Montgomery sub/add of residues in [0,p). */
mp_digit s_mp_ntt_mont_add(const ntt_mont *m, mp_digit a, mp_digit b)
{
   mp_digit d = a + b - m->p;
   return d + (m->p & ((mp_digit)0 - (d >> NTT_BASE_SIZE_M_ONE)));
}

mp_digit s_mp_ntt_mont_sub(const ntt_mont *m, mp_digit a, mp_digit b)
{
   mp_digit d = a - b;
   return d + (m->p & ((mp_digit)0 - (d >> NTT_BASE_SIZE_M_ONE)));
}

/* Convert x in [0,p) into Montgomery form (x * R mod p), and back. */
mp_digit s_mp_ntt_to_mont(const ntt_mont *m, mp_digit x)
{
   return s_mp_ntt_mont_mul(m, x, m->r2);
}

mp_digit s_mp_ntt_from_mont(const ntt_mont *m, mp_digit a)
{
   return s_mp_ntt_mont_redc(m, (mp_word)a);
}

/* TODO: we are down to three primes, no extra sort function necessary, do it inline */
/* Insertion sort by prime value. Both drivers select their primes largest
 * first and hand them to Garner ascending. */
void s_mp_ntt_sort_primes(ntt_prime_ent *ents, int n)
{
   int a, b;

   for (a = 1; a < n; a++) {
      ntt_prime_ent key = ents[a];
      b = a - 1;
      while (b >= 0 && ents[b].p > key.p) {
         ents[b + 1] = ents[b];
         b--;
      }
      ents[b + 1] = key;
   }
}

void s_mp_ntt_mont_init(ntt_mont *m, mp_digit p)
{
   mp_digit inv, r;
   m->p = p;
   inv = s_mp_ntt_invmod_d(p);
   m->ninv = (mp_digit)(0 - inv);                     /* -p^{-1} mod 2^64 */
   r = (mp_digit)(((mp_word)1 << NTT_BASE_SIZE) % p); /* R mod p = 2^64 mod p */
   m->r2 = (mp_digit)((mp_word)r * r % p);            /* R^2 mod p */
}

/*
 * Montgomery-form modular exponentiation: base^exp mod p, base in Montgomery
 * form, result in Montgomery form.
 */
mp_digit s_mp_ntt_mont_pow(const ntt_mont *m, mp_digit base, mp_digit exp)
{
   /* 1 in Montgomery form. No define because it is not a constant here */
   mp_digit result = s_mp_ntt_to_mont(m, 1);
   mp_digit b = base;
   while (exp) {
      if (exp & 1) result = s_mp_ntt_mont_mul(m, result, b);
      b = s_mp_ntt_mont_mul(m, b, b);
      exp >>= 1;
   }
   return result;
}

/*
 * Smallest 5-smooth number >= need.  The search enumerates 2^i * 3^j * 5^k
   directly; the bound 4*need is safe because consecutive 5-smooth numbers never
   differ by more than a factor 5/4 (for n > 8) and by much less once they
   are dense.
   The naive way (weights 1,1,1) is much slower on this machine.
 */
#ifdef USE_NAIVE
int s_mp_ntt_smooth_length_naive(int need)
{
   int64_t best = -1, p2, p3, p5;

   if (need <= 1) {
      return 1;
   }

   for (p2 = 1; p2 < (4LL * need); p2 *= 2) {
      for (p3 = p2; p3 < (4LL * need); p3 *= 3) {
         for (p5 = p3; p5 < (4LL * need); p5 *= 5) {
            if (p5 >= need && (best < 0 || p5 < best)) best = p5;
         }
      }
   }
   return (int)best;
}
#else
#include <limits.h>
/* Gives the different radix sizes a different weight according to runtime */
int s_mp_ntt_smooth_length_weighted(int n_min, int w2, int w3, int w5)
{
   int64_t next_pow2 = 1;
   int exp_pow2 = 0;
   int a, b, c, p5, p5_p3;
   int64_t n, current_n, size_limit, current_cost, cost;

   while (next_pow2 < n_min) {
      next_pow2 *= 2;
      exp_pow2++;
   }
   n = next_pow2;
   a = exp_pow2;
   b = 0;
   c = 0;
   cost = next_pow2 * (exp_pow2 * w2);

   size_limit = next_pow2 * 2;

   p5 = 1;
   for (c = 0; p5 <= size_limit; c++, p5 *= 5) {
      p5_p3 = p5;
      for (b = 0; p5_p3 <= size_limit; b++, p5_p3 *= 3) {
         /* Find the minimum 'a' to cross n_min */
         a = 0;
         current_n = p5_p3;
         while (current_n < n_min) {
            current_n *= 2;
            a++;
         }
         /* Just brute force */
         while (current_n <= size_limit) {
            current_cost = current_n * (a * w2 + b * w3 + c * w5);
            if (current_cost < cost) {
               cost = current_cost;
               n = current_n;
            }
            /* Pick the smaller in case of cost eqaulity */
            else if (current_cost == cost && current_n < n) {
               n = current_n;
            }
            current_n *= 2;
            a++;
         }
         /* Check for overflow */
         if (p5_p3 > (INT_MAX / 3)) {
            break;
         }
      }
      if (p5 > (INT_MAX / 5)) {
         break;
      }
   }
   return n;
}
#endif

int s_mp_ntt_smooth_length(int need)
{
   int length;
   /* TODO: was included for tuning. Either change to runtime or remove */
#ifdef USE_NAIVE
   length =  s_mp_ntt_smooth_length_naive(need);
#else
   length = s_mp_ntt_smooth_length_weighted(need, NTT_SMOOTH_W2, NTT_SMOOTH_W3, NTT_SMOOTH_W5);
#endif
   return length;
}


/*
   Radices of m in the order the stages apply them: all the 5s, then the 3s,
   then the 2s.  (Any order is correct; this one puts the long-radix work in the
   early stages, where the block is large and the twiddles are few or trivial,
   and leaves the lot of radix-2 stages for last, where a block is two adjacent
   elements.)  Returns the number of stages, or -1 if m is not 5-smooth or needs
   more room than the caller has.
 */
/* TODO: output to given variable and only return errors */
int s_mp_ntt_smooth_plan(int m, int *radix, int max)
{
   int n = 0, t = m;

   if (m < 2) return -1;
   while (t % 5 == 0) {
      if (n >= max) {
         return -1;
      }
      radix[n++] = 5;
      t /= 5;
   }
   while (t % 3 == 0) {
      if (n >= max) {
         return -1;
      }
      radix[n++] = 3;
      t /= 3;
   }
   while (t % 2 == 0) {
      if (n >= max) {
         return -1;
      }
      radix[n++] = 2;
      t /= 2;
   }
   return (t == 1) ? n : -1;
}


/*
 * Radix 5 needs more care than radix 2 or 3: the field has no imaginary part to
 * split the transform into.  With t = zeta^{m/5} and the pairs (y1,y4), (y2,y3),
 *
 *      t1 = y1+y4,  u1 = y1-y4,  t2 = y2+y3,  u2 = y2-y3,
 *      c1 = t + t^{-1},  c2 = t^2 + t^{-2},
 *      s1 = t - t^{-1},  s12 = t^2 - t^{-2} = s1*c1,
 *
 * the five outputs are
 *
 *      Z_0 = y0 + t1 + t2
 *      Z_1 = y0 + (t1*c1 + u1*s1 + t2*c2 + u2*s12)/2
 *      Z_4 = y0 + (t1*c1 - u1*s1 + t2*c2 - u2*s12)/2
 *      Z_2 = y0 + (t1*c2 + u1*s12 + t2*c1 - u2*s1)/2
 *      Z_3 = y0 + (t1*c2 - u1*s12 + t2*c1 + u2*s1)/2
 *
 * The 1/2 comes from writing each member of a pair through its sum and
 * difference; folding it into the constants leaves the butterfly with eight
 * multiplications per five elements - four distinct products per pair, each one
 * reused by a sum and a difference. The inverse uses the same constants with
 * s1 and s12 negated - those are the two that change sign under t -> t^{-1};
 * c1 and c2 are symmetric.
 */
mp_err s_mp_ntt_sroots_init(ntt_sroots *R, ntt_mont *mc, int m, mp_digit zeta)
{
   int i, ell, n;

   R->m = m;
   R->nstages = 0;
   R->w = NULL;
   R->has3 = ((m % 3) == 0);
   R->has5 = ((m % 5) == 0);
   R->c3 = R->c3i = 0;

   for (i = 0; i < 4; i++) {
      R->k5[i] = 0;
      R->k5i[i] = 0;
   }

   if (m < 2) {
      return MP_VAL;
   }
   n = s_mp_ntt_smooth_plan(m, R->radix, NTT_MR_MAXSTAGE);
   if (n <= 0) {
      return MP_VAL;
   }
   R->nstages = n;
   for (i = 0, ell = m; i < n; i++) {
      R->ell[i] = ell;
      ell /= R->radix[i];
   }

   R->w = (mp_digit *)MP_MALLOC((size_t)m * sizeof(mp_digit));
   if (!R->w) {
      return MP_MEM;
   }
   R->w[0] = s_mp_ntt_to_mont(mc, 1);
   for (i = 1; i < m; i++) {
      R->w[i] = s_mp_ntt_mont_mul(mc, R->w[i - 1], zeta);
   }
   if (R->has3) {
      R->c3  = R->w[m / 3];            /* primitive cube root */
      R->c3i = R->w[m - m / 3];        /* its inverse */
   }
   if (R->has5) {
      mp_digit t   = R->w[m / 5];      /* primitive fifth root */
      mp_digit ti  = R->w[m - m / 5];  /* its inverse */
      mp_digit h   = s_mp_ntt_to_mont(mc, (mc->p + 1) >> 1);   /* 2^{-1} mod p */
      mp_digit c1  = s_mp_ntt_mont_add(mc, t, ti);
      /* TODO: is that readable or should distinct temporary variables be used? */
      mp_digit c2  = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_mul(mc, t, t), s_mp_ntt_mont_mul(mc, ti, ti));
      mp_digit s1  = s_mp_ntt_mont_sub(mc, t, ti);
      mp_digit s12 = s_mp_ntt_mont_mul(mc, s1, c1);

      R->k5[0] = s_mp_ntt_mont_mul(mc, c1,  h);
      R->k5[1] = s_mp_ntt_mont_mul(mc, c2,  h);
      R->k5[2] = s_mp_ntt_mont_mul(mc, s1,  h);
      R->k5[3] = s_mp_ntt_mont_mul(mc, s12, h);
      R->k5i[0] = R->k5[0];
      R->k5i[1] = R->k5[1];
      R->k5i[2] = s_mp_ntt_mont_sub(mc, 0, R->k5[2]);   /* -h*s1  */
      R->k5i[3] = s_mp_ntt_mont_sub(mc, 0, R->k5[3]);   /* -h*s12 */
   }
   return MP_OKAY;
}


void s_mp_ntt_sroots_free(ntt_sroots *R)
{
   if (!R) {
      return;
   }
   MP_FREE(R->w, (size_t)R->m * sizeof(mp_digit));
   R->w = NULL;
   R->nstages = 0;
}

/*
 * Forward: transform first and do the twiddles after.
 * Inverse: undoes the twiddles first (reading the same table at m - e,
 * since omega^{-e} = omega^{m-e}) and inverse transform after.
 * The inverse is unscaled, because s_mp_ntt_smooth_inv divides by m once at
 * the end. Slot q carries index q in both directions.
 */
void s_mp_ntt_fwd2(const ntt_mont *mc, mp_digit *b, int L, int step, const ntt_sroots *R)
{
   const mp_digit *w = R->w;
   int m = R->m, e = 0, j;

   for (j = 0; j < L; j++) {
      mp_digit y0 = b[j], y1 = b[j + L];
      b[j]     = s_mp_ntt_mont_add(mc, y0, y1);                    /* q = 0, W^0 = 1 */
      b[j + L] = s_mp_ntt_mont_mul(mc, s_mp_ntt_mont_sub(mc, y0, y1), w[e]);
      e += step;
      if (e >= m) {
         e -= m;
      }
   }
}

void s_mp_ntt_fwd3(const ntt_mont *mc, mp_digit *b, int L, int step, const ntt_sroots *R)
{
   const mp_digit *w = R->w;
   int m = R->m, j;
   int s1 = step, s2 = 2 * step;       /* step = m/ell <= m/3, so 2*step < m */
   int e1 = 0, e2 = 0;
   mp_digit t = R->c3;

   for (j = 0; j < L; j++) {
      mp_digit y0 = b[j], y1 = b[j + L], y2 = b[j + 2 * L];
      mp_digit v  = s_mp_ntt_mont_mul(mc, s_mp_ntt_mont_sub(mc, y1, y2), t);          /* t*(y1-y2) */
      mp_digit z0 = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, y1), y2);
      mp_digit z1 = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_sub(mc, y0, y2), v);
      mp_digit z2 = s_mp_ntt_mont_sub(mc, s_mp_ntt_mont_sub(mc, y0, y1), v);

      b[j]         = z0;
      b[j + L]     = s_mp_ntt_mont_mul(mc, z1, w[e1]);
      b[j + 2 * L] = s_mp_ntt_mont_mul(mc, z2, w[e2]);
      e1 += s1;
      if (e1 >= m) {
         e1 -= m;
      }
      e2 += s2;
      if (e2 >= m) {
         e2 -= m;
      }
   }
}

void s_mp_ntt_fwd5(const ntt_mont *mc, mp_digit *b, int L, int step, const ntt_sroots *R)
{
   const mp_digit *w = R->w;
   int m = R->m, j;
   int s1 = step, s2 = 2 * step, s3 = 3 * step, s4 = 4 * step;   /* all < m */
   int e1 = 0, e2 = 0, e3 = 0, e4 = 0;
   mp_digit k0 = R->k5[0], k1 = R->k5[1], k2 = R->k5[2], k3 = R->k5[3];
   mp_digit y0, y1, y2, y3, y4, t1, t2, u1, u2, t1c1, t2c2, t1c2, t2c1, u1s1, u2s12, u1s12, u2s1, z1, z2, z3, z4;

   for (j = 0; j < L; j++) {
      y0 = b[j];
      y1 = b[j + L];
      y2 = b[j + 2 * L];
      y3 = b[j + 3 * L];
      y4 = b[j + 4 * L];
      t1 = s_mp_ntt_mont_add(mc, y1, y4);
      u1 = s_mp_ntt_mont_sub(mc, y1, y4);
      t2 = s_mp_ntt_mont_add(mc, y2, y3);
      u2 = s_mp_ntt_mont_sub(mc, y2, y3);
      t1c1 = s_mp_ntt_mont_mul(mc, t1, k0);
      u1s1  = s_mp_ntt_mont_mul(mc, u1, k2);
      t2c2 = s_mp_ntt_mont_mul(mc, t2, k1);
      u2s12 = s_mp_ntt_mont_mul(mc, u2, k3);
      t1c2 = s_mp_ntt_mont_mul(mc, t1, k1);
      u1s12 = s_mp_ntt_mont_mul(mc, u1, k3);
      t2c1 = s_mp_ntt_mont_mul(mc, t2, k0);
      u2s1  = s_mp_ntt_mont_mul(mc, u2, k2);
      z1 = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, s_mp_ntt_mont_add(mc, t1c1, u1s1)),  s_mp_ntt_mont_add(mc, t2c2,
                             u2s12));
      z2 = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, s_mp_ntt_mont_add(mc, t1c2, u1s12)), s_mp_ntt_mont_sub(mc, t2c1,
                             u2s1));
      z3 = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, s_mp_ntt_mont_sub(mc, t1c2, u1s12)), s_mp_ntt_mont_add(mc, t2c1,
                             u2s1));
      z4 = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, s_mp_ntt_mont_sub(mc, t1c1, u1s1)),  s_mp_ntt_mont_sub(mc, t2c2,
                             u2s12));

      b[j]         = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, t1), t2);
      b[j + L]     = s_mp_ntt_mont_mul(mc, z1, w[e1]);
      b[j + 2 * L] = s_mp_ntt_mont_mul(mc, z2, w[e2]);
      b[j + 3 * L] = s_mp_ntt_mont_mul(mc, z3, w[e3]);
      b[j + 4 * L] = s_mp_ntt_mont_mul(mc, z4, w[e4]);
      e1 += s1;
      if (e1 >= m) {
         e1 -= m;
      }
      e2 += s2;
      if (e2 >= m) {
         e2 -= m;
      }
      e3 += s3;
      if (e3 >= m) {
         e3 -= m;
      }
      e4 += s4;
      if (e4 >= m) {
         e4 -= m;
      }
   }
}

/* INTTs */
void s_mp_ntt_inv2(const ntt_mont *mc, mp_digit *b, int L, int step, const ntt_sroots *R)
{
   const mp_digit *w = R->w;
   int m = R->m, e = 0, j;

   for (j = 0; j < L; j++) {
      mp_digit y0 = b[j];
      mp_digit y1 = s_mp_ntt_mont_mul(mc, b[j + L], w[e ? m - e : 0]);
      b[j]     = s_mp_ntt_mont_add(mc, y0, y1);
      b[j + L] = s_mp_ntt_mont_sub(mc, y0, y1);
      e += step;
      if (e >= m) {
         e -= m;
      }
   }
}

void s_mp_ntt_inv3(const ntt_mont *mc, mp_digit *b, int L, int step, const ntt_sroots *R)
{
   const mp_digit *w = R->w;
   int m = R->m, j;
   int s1 = step, s2 = 2 * step;
   int e1 = 0, e2 = 0;
   mp_digit t = R->c3i;

   for (j = 0; j < L; j++) {
      mp_digit y0 = b[j];
      mp_digit y1 = s_mp_ntt_mont_mul(mc, b[j + L],     w[e1 ? m - e1 : 0]);
      mp_digit y2 = s_mp_ntt_mont_mul(mc, b[j + 2 * L], w[e2 ? m - e2 : 0]);
      mp_digit v  = s_mp_ntt_mont_mul(mc, s_mp_ntt_mont_sub(mc, y1, y2), t);

      b[j]         = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, y1), y2);
      b[j + L]     = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_sub(mc, y0, y2), v);
      b[j + 2 * L] = s_mp_ntt_mont_sub(mc, s_mp_ntt_mont_sub(mc, y0, y1), v);
      e1 += s1;
      if (e1 >= m) {
         e1 -= m;
      }
      e2 += s2;
      if (e2 >= m) {
         e2 -= m;
      }
   }
}

void s_mp_ntt_inv5(const ntt_mont *mc, mp_digit *b, int L, int step, const ntt_sroots *R)
{
   const mp_digit *w = R->w;
   int m = R->m, j;
   int s1 = step, s2 = 2 * step, s3 = 3 * step, s4 = 4 * step;
   int e1 = 0, e2 = 0, e3 = 0, e4 = 0;
   mp_digit k0 = R->k5i[0], k1 = R->k5i[1], k2 = R->k5i[2], k3 = R->k5i[3];
   mp_digit y0, y1, y2, y3, y4, t1, t2, u1, u2, t1c1, t2c2, t1c2, t2c1, u1s1, u2s12, u1s12, u2s1;


   for (j = 0; j < L; j++) {
      y0 = b[j];
      y1 = s_mp_ntt_mont_mul(mc, b[j + L],     w[e1 ? m - e1 : 0]);
      y2 = s_mp_ntt_mont_mul(mc, b[j + 2 * L], w[e2 ? m - e2 : 0]);
      y3 = s_mp_ntt_mont_mul(mc, b[j + 3 * L], w[e3 ? m - e3 : 0]);
      y4 = s_mp_ntt_mont_mul(mc, b[j + 4 * L], w[e4 ? m - e4 : 0]);
      t1 = s_mp_ntt_mont_add(mc, y1, y4);
      u1 = s_mp_ntt_mont_sub(mc, y1, y4);
      t2 = s_mp_ntt_mont_add(mc, y2, y3);
      u2 = s_mp_ntt_mont_sub(mc, y2, y3);
      t1c1 = s_mp_ntt_mont_mul(mc, t1, k0);
      u1s1  = s_mp_ntt_mont_mul(mc, u1, k2);
      t2c2 = s_mp_ntt_mont_mul(mc, t2, k1);
      u2s12 = s_mp_ntt_mont_mul(mc, u2, k3);
      t1c2 = s_mp_ntt_mont_mul(mc, t1, k1);
      u1s12 = s_mp_ntt_mont_mul(mc, u1, k3);
      t2c1 = s_mp_ntt_mont_mul(mc, t2, k0);
      u2s1  = s_mp_ntt_mont_mul(mc, u2, k2);

      b[j]         = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, t1), t2);
      b[j + L]     = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, s_mp_ntt_mont_add(mc, t1c1, u1s1)),
                                       s_mp_ntt_mont_add(mc, t2c2, u2s12));
      b[j + 2 * L] = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, s_mp_ntt_mont_add(mc, t1c2, u1s12)),
                                       s_mp_ntt_mont_sub(mc, t2c1, u2s1));
      b[j + 3 * L] = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, s_mp_ntt_mont_sub(mc, t1c2, u1s12)),
                                       s_mp_ntt_mont_add(mc, t2c1, u2s1));
      b[j + 4 * L] = s_mp_ntt_mont_add(mc, s_mp_ntt_mont_add(mc, y0, s_mp_ntt_mont_sub(mc, t1c1, u1s1)),
                                       s_mp_ntt_mont_sub(mc, t2c2, u2s12));
      e1 += s1;
      if (e1 >= m) {
         e1 -= m;
      }
      e2 += s2;
      if (e2 >= m) {
         e2 -= m;
      }
      e3 += s3;
      if (e3 >= m) {
         e3 -= m;
      }
      e4 += s4;
      if (e4 >= m) {
         e4 -= m;
      }
   }
}

void s_mp_ntt_smooth_fwd(ntt_mont *mc, mp_digit *a, int m, const ntt_sroots *R)
{
   int k, start, r, ell, L, step;

   for (k = 0; k < R->nstages; k++) {
      r = R->radix[k];
      ell = R->ell[k];
      L = ell / r;
      step = m / ell;
      for (start = 0; start < m; start += ell) {
         mp_digit *b = a + start;
         switch (r) {
         case 5:
            s_mp_ntt_fwd5(mc, b, L, step, R);
            break;
         case 3:
            s_mp_ntt_fwd3(mc, b, L, step, R);
            break;
         default:
            s_mp_ntt_fwd2(mc, b, L, step, R);
            break;
         }
      }
   }
}

void s_mp_ntt_smooth_inv(ntt_mont *mc, mp_digit *a, int m, const ntt_sroots *R)
{
   int k, start, i;
   mp_digit ninv;

   for (k = R->nstages - 1; k >= 0; k--) {
      int r = R->radix[k], ell = R->ell[k], L = ell / r, step = m / ell;
      for (start = 0; start < m; start += ell) {
         mp_digit *b = a + start;
         switch (r) {
         case 5:
            s_mp_ntt_inv5(mc, b, L, step, R);
            break;
         case 3:
            s_mp_ntt_inv3(mc, b, L, step, R);
            break;
         default:
            s_mp_ntt_inv2(mc, b, L, step, R);
            break;
         }
      }
   }
   ninv = s_mp_ntt_mont_pow(mc, s_mp_ntt_to_mont(mc, (mp_digit)m), mc->p - 2);   /* 1/m, Montgomery */
   for (i = 0; i < m; i++) {
      a[i] = s_mp_ntt_mont_mul(mc, a[i], ninv);
   }
}

/* TODO: Is it somewhere else, too? Look for it. */
int s_mp_ntt_ilog2_ceil(int x)
{
   int s = 0;
   while ((1 << s) < x) {
      s++;
   }
   return s;
}

/*
 * One prime's contribution: forward-transform both inputs, multiply pointwise,
 * inverse-transform, and return the convolution coefficients modulo that prime.
 *
 * The transform runs at m = ent->m (5-smooth, >= la + lb - 1), so its cyclic
 * convolution of length m agrees with the linear one on the first
 * ncoef = la + lb - 1 coefficients and is zero beyond - which is why only those
 * are read back.
 */
mp_err s_mp_ntt_convolve_single(const mp_int *a, const mp_int *b, const ntt_prime_ent *ent, mp_digit *ra,
                                mp_digit *rb, mp_digit *out, int ncoef)
{
   int m = ent->m, i;
   ntt_mont mc;
   ntt_sroots R;
   mp_digit zeta;
   mp_err err = MP_OKAY;

   s_mp_ntt_mont_init(&mc, ent->p);
   zeta = s_mp_ntt_to_mont(&mc, ent->root);        /* primitive m-th root (Montgomery) */

   if ((err = s_mp_ntt_sroots_init(&R, &mc, m, zeta)) != MP_OKAY) return err;

   for (i = 0; i < m; i++) {
      ra[i] = 0;
      rb[i] = 0;
   }
   for (i = 0; i < a->used; i++) {
      ra[i] = s_mp_ntt_to_mont(&mc, a->dp[i]);
   }
   for (i = 0; i < b->used; i++) {
      rb[i] = s_mp_ntt_to_mont(&mc, b->dp[i]);
   }

   s_mp_ntt_smooth_fwd(&mc, ra, m, &R);
   s_mp_ntt_smooth_fwd(&mc, rb, m, &R);

   for (i = 0; i < m; i++) {
      ra[i] = s_mp_ntt_mont_mul(&mc, ra[i], rb[i]);
   }

   s_mp_ntt_smooth_inv(&mc, ra, m, &R);            /* includes the 1/m */

   for (i = 0; i < ncoef; i++) {
      out[i] = s_mp_ntt_from_mont(&mc, ra[i]);
   }

   s_mp_ntt_sroots_free(&R);
   return err;
}

mp_err s_mp_ntt_autoconvolve_single(const mp_int *a, const ntt_prime_ent *ent, mp_digit *ra, mp_digit *out, int ncoef)
{
   int m = ent->m, i;
   ntt_mont mc;
   ntt_sroots R;
   mp_digit zeta;
   mp_err err = MP_OKAY;

   s_mp_ntt_mont_init(&mc, ent->p);
   zeta = s_mp_ntt_to_mont(&mc, ent->root);        /* primitive m-th root (Montgomery) */

   if ((err = s_mp_ntt_sroots_init(&R, &mc, m, zeta)) != MP_OKAY) return err;

   for (i = 0; i < m; i++) {
      ra[i] = 0;
   }
   for (i = 0; i < a->used; i++) {
      ra[i] = s_mp_ntt_to_mont(&mc, a->dp[i]);
   }

   s_mp_ntt_smooth_fwd(&mc, ra, m, &R);

   for (i = 0; i < m; i++) {
      ra[i] = s_mp_ntt_mont_mul(&mc, ra[i], ra[i]);
   }

   s_mp_ntt_smooth_inv(&mc, ra, m, &R);            /* includes the 1/m */

   for (i = 0; i < ncoef; i++) {
      out[i] = s_mp_ntt_from_mont(&mc, ra[i]);
   }

   s_mp_ntt_sroots_free(&R);
   return err;
}



mp_err s_mp_ntt_mul_smooth(const mp_int *a, const mp_int *b, mp_int *c)
{
   int ncoef, m, minl, bound_bits, nprimes, pi, i;
   ntt_prime_ent *ents = NULL;
   ntt_crt_ctx   *cx = NULL;
   mp_digit *primes = NULL, **residues = NULL, *ra = NULL, *rb = NULL;
   mp_int P, pmp;
   mp_err err = MP_OKAY;

   if (mp_iszero(a) || mp_iszero(b)) {
      return mp_mul(a, b, c);
   }
   if (a->used < 8 || b->used < 8) {
      return mp_mul(a, b, c);
   }

   ncoef = a->used + b->used - 1; /* coefficients that carry the product */
   if ((ncoef > (1 << NTT_QUARTER_COEFF)) || (ncoef < 0)) {
      return MP_VAL;
   }
   m = s_mp_ntt_smooth_length(ncoef); /* smallest 5-smooth transform length */
   if (m < ncoef) {
      return MP_VAL;
   }

   /* Convolution bound: each output coefficient < min(la,lb) * 2^120. */
   minl = (a->used < b->used) ? a->used : b->used;
   bound_bits = (2 * NTT_LIMB_SIZE) + s_mp_ntt_ilog2_ceil(minl) + 1;
   if (((bound_bits + NTT_LIMB_SIZE_M_ONE) / NTT_LIMB_SIZE) > NTT_MAXW) {
      return MP_VAL;
   }

   ents = (ntt_prime_ent *)MP_MALLOC((size_t)NTT_MAXP * sizeof(ntt_prime_ent));
   if (ents == NULL) {
      return MP_MEM;
   }

   if ((err = mp_init_multi(&P, &pmp, NULL)) != MP_OKAY) {
      MP_FREE(ents, (size_t)NTT_MAXP * sizeof(ntt_prime_ent));
      return err;
   }
   /* TODO */
   mp_set(&P, 1);

   /* Find primes until their product > 2^bound_bits. should result in three primes */
   nprimes = 0;
   while (mp_count_bits(&P) <= bound_bits) {
      if (nprimes >= NTT_MAXP) {
         err = MP_VAL;
         goto LTM_ERR;
      }
      if (!s_mp_ntt_find_smooth_prime(&ents[nprimes], m, nprimes, 0)) {
         err = MP_VAL;
         goto LTM_ERR;
      }
      mp_set(&pmp, ents[nprimes].p);
      if ((err = mp_mul(&P, &pmp, &P)) != MP_OKAY) goto LTM_ERR;
      nprimes++;
   }


   s_mp_ntt_sort_primes(ents, nprimes);

   ra = (mp_digit *)MP_MALLOC((size_t)m * sizeof(mp_digit));
   if (ra == NULL) {
      err = MP_MEM;
      goto LTM_ERR;
   }
   rb = (mp_digit *)MP_MALLOC((size_t)m * sizeof(mp_digit));
   if (rb == NULL) {
      err = MP_MEM;
      goto LTM_ERR;
   }
   residues = (mp_digit **)MP_CALLOC((size_t)nprimes, sizeof(mp_digit *));
   if (residues == NULL) {
      err = MP_MEM;
      goto LTM_ERR;
   }
   primes = (mp_digit *)MP_MALLOC((size_t)nprimes * sizeof(mp_digit));
   if (primes == NULL) {
      err = MP_MEM;
      goto LTM_ERR;
   }


   /* convolve over each prime and keep the residues */
   for (pi = 0; pi < nprimes; pi++) {
      residues[pi] = (mp_digit *)MP_MALLOC((size_t)ncoef * sizeof(mp_digit));
      if (residues[pi] == NULL) {
         err = MP_MEM;
         goto LTM_ERR;
      }
      if ((err = s_mp_ntt_convolve_single(a, b, &ents[pi], ra, rb, residues[pi], ncoef)) != MP_OKAY)        goto LTM_ERR;
      primes[pi] = ents[pi].p;
   }

   /* CRT: residues -> binary, then one linear carry sweep */
   cx = (ntt_crt_ctx *)MP_MALLOC(sizeof(ntt_crt_ctx));
   if (cx == NULL) {
      err = MP_MEM;
      goto LTM_ERR;
   }
   s_mp_ntt_crt_init(cx, primes, nprimes, bound_bits);
   err = s_mp_ntt_crt_assemble(cx, residues, nprimes, ncoef, c);

LTM_ERR:
   mp_clear_multi(&P, &pmp, NULL);
   if (residues) {
      for (i = 0; i < nprimes; i++) {
         MP_FREE(residues[i], (size_t)ncoef * sizeof(mp_digit));
      }
      MP_FREE(residues, (size_t)nprimes * sizeof(mp_digit *));
   }
   MP_FREE(primes, (size_t)nprimes * sizeof(mp_digit));
   MP_FREE(ra, (size_t)m * sizeof(mp_digit));
   MP_FREE(rb, (size_t)m * sizeof(mp_digit));
   MP_FREE(cx, sizeof(ntt_crt_ctx));
   MP_FREE(ents, (size_t)NTT_MAXP * sizeof(ntt_prime_ent));
   return err;
}

mp_err s_mp_ntt_sqr_smooth(const mp_int *a, mp_int *c)
{
   mp_err err = MP_OKAY;
   /* TODO: nprimes is three, defined in NTT_MAXP, rewrite accordingly */
   int ncoef, m, minl, bound_bits, nprimes, pi, i;
   ntt_prime_ent *ents = NULL;
   ntt_crt_ctx   *cx = NULL;
   mp_digit *primes = NULL, **residues = NULL, *ra = NULL;
   mp_int P, pmp;


   if (mp_iszero(a)) {
      return mp_sqr(a,c);
   }
   if (a->used < 8) {
      return mp_sqr(a, c);
   }

   ncoef = a->used + a->used - 1; /* coefficients that carry the product */
   if (ncoef > (1 << NTT_QUARTER_COEFF) || ncoef < 0) {
      return MP_VAL;
   }
   m = s_mp_ntt_smooth_length(ncoef); /* smallest 5-smooth transform length */
   if (m < ncoef) {
      return MP_VAL;
   }

   /* Convolution bound: each output coefficient < min(la,lb) * 2^120. */
   minl = a->used;
   bound_bits = NTT_LARGE_SIZE + s_mp_ntt_ilog2_ceil(minl) + 1;
   if ((bound_bits + NTT_LIMB_SIZE_M_ONE) / NTT_LIMB_SIZE > NTT_MAXW) {
      return MP_VAL;
   }

   ents = (ntt_prime_ent *)MP_MALLOC((size_t)NTT_MAXP * sizeof(ntt_prime_ent));
   if (ents == NULL) {
      return MP_MEM;
   }

   if ((err = mp_init_multi(&P, &pmp, NULL)) != MP_OKAY) {
      MP_FREE(ents, (size_t)NTT_MAXP * sizeof(ntt_prime_ent));
      return err;
   }
   /* TODO */
   mp_set(&P, 1);

   /* Find primes until their product > 2^bound_bits. should result in three primes */
   nprimes = 0;
   while (mp_count_bits(&P) <= bound_bits) {
      if (nprimes >= NTT_MAXP) {
         err = MP_VAL;
         goto LTM_ERR;
      }
      if (!s_mp_ntt_find_smooth_prime(&ents[nprimes], m, nprimes, 0)) {
         err = MP_VAL;
         goto LTM_ERR;
      }
#ifdef MP_PRINT_NTT_PRIMES
      fprintf(stderr,"Prime nprimes %d %"PRIu64"\n", nprimes, ents[nprimes].p);
#endif
      mp_set(&pmp, ents[nprimes].p);
      if ((err = mp_mul(&P, &pmp, &P)) != MP_OKAY) goto LTM_ERR;
      nprimes++;
   }

   s_mp_ntt_sort_primes(ents, nprimes);

   ra = (mp_digit *)MP_MALLOC((size_t)m * sizeof(mp_digit));
   if (ra == NULL) {
      err = MP_MEM;
      goto LTM_ERR;
   }
   residues = (mp_digit **)MP_CALLOC((size_t)nprimes, sizeof(mp_digit *));
   if (residues == NULL) {
      err = MP_MEM;
      goto LTM_ERR;
   }
   primes = (mp_digit *)MP_MALLOC((size_t)nprimes * sizeof(mp_digit));
   if (primes == NULL) {
      err = MP_MEM;
      goto LTM_ERR;
   }


   /* convolve over each prime and keep the residues */
   for (pi = 0; pi < nprimes; pi++) {
      residues[pi] = (mp_digit *)MP_MALLOC((size_t)ncoef * sizeof(mp_digit));
      if (residues[pi] == NULL) {
         err = MP_MEM;
         goto LTM_ERR;
      }
      if ((err = s_mp_ntt_autoconvolve_single(a, &ents[pi], ra, residues[pi], ncoef)) != MP_OKAY)        goto LTM_ERR;
      primes[pi] = ents[pi].p;
   }

   /* CRT: residues -> binary, then one linear carry sweep */
   cx = (ntt_crt_ctx *)MP_MALLOC(sizeof(ntt_crt_ctx));
   if (cx == NULL) {
      err = MP_MEM;
      goto LTM_ERR;
   }
   s_mp_ntt_crt_init(cx, primes, nprimes, bound_bits);
   err = s_mp_ntt_crt_assemble(cx, residues, nprimes, ncoef, c);

LTM_ERR:
   mp_clear_multi(&P, &pmp, NULL);
   if (residues) {
      for (i = 0; i < nprimes; i++) {
         MP_FREE(residues[i], (size_t)ncoef * sizeof(mp_digit));
      }
      MP_FREE(residues, (size_t)nprimes * sizeof(mp_digit *));
   }
   MP_FREE(primes, (size_t)nprimes * sizeof(mp_digit));
   MP_FREE(ra, (size_t)m * sizeof(mp_digit));
   MP_FREE(cx, sizeof(ntt_crt_ctx));
   MP_FREE(ents, (size_t)NTT_MAXP * sizeof(ntt_prime_ent));
   return err;
}


/*
 * Search for the `index`-th largest usable prime of a 5-smooth transform length
 * m, i.e. p = k*m + 1 < 2^62 with k walked downwards (index = 0 is the largest),
 * requiring p > min_p.
 *
 * Largest rather than smallest for the same reason as ntt_find_large_prime:
 * three ~62-bit primes cover the whole convolution bound at every supported
 * size, so the prime count never grows with the input.
 *
 * The root search differs.  For the truncated path a quadratic non-residue gives
 * a 2^s-th root for free (Euler's criterion), but m has 3s and 5s in it and no
 * single test covers those.  So: w = a^((p-1)/m) has order dividing m, and it is
 * exactly m as soon as w^(m/r) != 1 for every prime r dividing m - because every
 * proper divisor of m divides one of those m/r.
 */
bool s_mp_ntt_find_smooth_prime(ntt_prime_ent *ent, int m, int index, mp_digit min_p)
{
   mp_digit mm, kmax, root;
   int seen = 0;
   mp_digit k, a, p, w;

   if (m < 2) {
      return false;
   }
   mm   = (mp_digit)m;
   kmax = (((mp_digit)1 << NTT_PRIME_MAX_SIZE) - 1) / mm;

   for (k = kmax; k >= 1; k--) {
      p = k * mm + 1;
      if (p >= ((mp_digit)1 << NTT_PRIME_MAX_SIZE)) {
         continue;
      }
      if (!s_mp_ntt_is_prime(p)) {
         continue;
      }
      if (p <= min_p) {
         return false;  /* cannot cover the bound */
      }
      if (seen++ != index) {
         continue;  /* want the index-th largest */
      }

      root = 0;
      for (a = 2; (a < 1000) && (root == 0); a++) {
         w = s_mp_small_powmod(a, k, p);             /* order divides m */
         if (w == 1) {
            continue;
         }
         if (((m % 2) == 0) && s_mp_small_powmod(w, (mp_digit)m / 2, p) == 1) {
            continue;
         }
         if (((m % 3) == 0) && s_mp_small_powmod(w, (mp_digit)m / 3, p) == 1) {
            continue;
         }
         if (((m % 5) == 0) && s_mp_small_powmod(w, (mp_digit)m / 5, p) == 1) {
            continue;
         }
         root = w;
      }
      if (root == 0) {
         return false;  /* no root of full order */
      }
      ent->p = p;
      ent->root = root;
      ent->s = 0;
      ent->m = m;
      return true;
   }
   return false;
}



void s_mp_ntt_crt_init(ntt_crt_ctx *cx, const mp_digit *primes, int t, int bound_bits)
{
   int i, j;

   cx->t = t;
   cx->width = (bound_bits + NTT_LIMB_SIZE_M_ONE) / NTT_LIMB_SIZE;
   if (cx->width < 1) cx->width = 1;

   for (i = 0; i < t; i++) {
      cx->p[i] = primes[i];
      s_mp_ntt_mont_init(&cx->mc[i], primes[i]);
   }
   for (i = 0; i < t; i++) {
      for (j = i + 1; j < t; j++) {
         mp_digit pj = cx->p[j];
         mp_digit inv = s_mp_small_powmod(cx->p[i] % pj, pj - 2, pj);  /* p_i^{-1} mod p_j */
         cx->minv[i][j] = s_mp_ntt_to_mont(&cx->mc[j], inv);
      }
   }
}

/* Recover one coefficient (from its residue modulo every prime) as a short
 * array of 60-bit limbs.
 *
 * minv[i][j] = p_i^{-1} * 2^64 mod p_j, for i < j.  Propagating the digit d_i to
 * the higher residues divides by p_i, so the update is
 *     x_j <- (x_j - d_i) * p_i^{-1} mod p_j.
 * (M_i^{-1} = (p_0 ... p_{i-1})^{-1}, not p_i^{-1}, appears only when the partial
 * value rather than the digit is propagated.)  Keeping p_i^{-1} in Montgomery
 * form makes that update one multiplication with no division.
 */
void s_mp_ntt_crt_coeff(ntt_crt_ctx *cx, const mp_digit *r, mp_digit *limbs)
{
   int t = cx->t, i, j, w;
   mp_digit *x = cx->x;
   mp_digit *d = cx->d;

   for (i = 0; i < t; i++) x[i] = r[i];        /* each < p_i */

   for (i = 0; i < t; i++) {
      mp_digit di = x[i];                          /* mixed-radix digit */
      d[i] = di;
      for (j = i + 1; j < t; j++) {
         mp_digit pj = cx->p[j];
         mp_digit diff = (x[j] >= di) ? (x[j] - di) : (x[j] + pj - di);
         x[j] = s_mp_ntt_mont_mul(&cx->mc[j], diff, cx->minv[i][j]);
      }
   }

   /* Horner: c = d_0 + p_0 (d_1 + p_1 (d_2 + ... )) */
   for (w = 0; w < NTT_MAXW; w++) limbs[w] = 0;
   limbs[0] = d[t - 1];
   for (i = t - 2; i >= 0; i--) {
      mp_word carry = (mp_word)d[i];
      for (w = 0; w < NTT_MAXW; w++) {
         mp_word v = (mp_word)limbs[w] * cx->p[i] + carry;
         limbs[w] = (mp_digit)(v & NTT_MASK);
         carry = v >> NTT_LIMB_SIZE;
      }
   }
}

/* Residues -> the product.  residues[i] holds the first ncoef convolution
 * coefficients modulo primes[i]; the result is written into c. */
mp_err s_mp_ntt_crt_assemble(ntt_crt_ctx *cx, mp_digit **residues, int nprimes,
                             int ncoef, mp_int *c)
{
   mp_err err = MP_OKAY;
   mp_digit *r, *acc, limbs[NTT_MAXW];
   mp_digit carry;
   int k, i, w, n;
   mp_word v;

   r = (mp_digit *)MP_MALLOC((size_t)nprimes * sizeof(mp_digit));
   if (r == NULL) {
      return MP_MEM;
   }
   /* lazy accumulator: one 60-bit slot per limb position; a slot receives at
    * most one limb from each of the <= NTT_MAXW overlapping coefficients, so it
    * stays below NTT_MAXW * 2^60 and cannot overflow */
   acc = (mp_digit *)MP_CALLOC((size_t)ncoef + NTT_MAXW, sizeof(mp_digit));
   if (acc == NULL) {
      MP_FREE(r, (size_t)nprimes * sizeof(mp_digit));
      return MP_MEM;
   }

   for (k = 0; k < ncoef; k++) {
      for (i = 0; i < nprimes; i++) {
         r[i] = residues[i][k];
      }
      s_mp_ntt_crt_coeff(cx, r, limbs);
      for (w = 0; w < cx->width; w++) {
         acc[k + w] += limbs[w];
      }
   }

   n = ncoef + NTT_MAXW;
   if ((err = mp_grow(c, n)) != MP_OKAY)        goto LTM_ERR;
   if (err == MP_OKAY) {
      carry = 0;
      for (k = 0; k < n; k++) {
         v = (mp_word)acc[k] + carry;
         c->dp[k] = (mp_digit)(v & NTT_MASK);
         carry = (mp_digit)(v >> NTT_LIMB_SIZE);
      }
      c->used = n;
      c->sign = MP_ZPOS;
      mp_clamp(c);
   }
LTM_ERR:
   MP_FREE(r, (size_t)nprimes * sizeof(mp_digit));
   MP_FREE(acc, ((size_t)ncoef + NTT_MAXW) * sizeof(mp_digit));
   return err;
}
#else
#  ifdef _MSC_VER
#    pragma message("NTT for MP_32BIT and MP_64BIT only.")
#  else
#    warning "NTT for MP_32BIT and MP_64BIT only."
#  endif
#endif

#endif
