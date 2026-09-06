#include "tommath_private.h"
#ifdef MP_FACTORIAL_DIVISORS_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

/*
   Legendre's formula (also called de Polignac's formula).

      E_{p}(n!)=\sum _{k=1}^{\infty }\left\lfloor \frac{n}{p^{k}}\right\rfloor

   It calculates the exponent of a prime $p$ in the prime factorization of
   $n!$ (n factorial). In number theory, this is written as $v_p(n!)$, which
   represents the highest power of $p$ that divides $n!$.
*/
mp_err mp_factorial_divisors(const mp_digit n, const mp_digit p, mp_digit *d)
{
   mp_digit q, m;
   m = 0u;
   if ((p < 2) || (n == 0)) {
      *d = 0u;
      return MP_VAL;
   }
   if (p > n) {
      *d = 0u;
      return MP_OKAY;
   }
   if (p > (n / 2u)) {
      *d = 1u;
      return MP_OKAY;
   }
   q = n;
   while (q >= p) {
      q = q / p;
      m += q;
   }
   *d = m;
   return MP_OKAY;
}

#endif
