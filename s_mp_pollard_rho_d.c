#include "tommath_private.h"
#ifdef S_MP_POLLARD_RHO_D_C
/* LibTomMath, multiple-precision integer library -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */



mp_digit s_mp_pollard_rho_d(mp_digit n, int limit)
{
   mp_digit startvalue = 2;
   mp_digit x, y, t;
   mp_digit d = 1;
   int i = 0;

   /* srand(0xdeadbeef); */

   x = y = startvalue;
   while (d == 1) {
      x = s_mp_sqrincrmod_d(x, n);
      t = s_mp_sqrincrmod_d(y, n);
      y = s_mp_sqrincrmod_d(t, n);
      d = (y > x) ? y - x : x - y;
      d = s_mp_gcd_d(d,n);
      if (d == n) {
         startvalue++;
         x = y = startvalue;
      }
      if ((i % 10) == 0) {
         if (i >= limit) {
            fprintf(stderr,"s_mp_pollard_rho reached limit\n");
            return 0u;
         }
         /* startvalue = (mp_digit)rand(); */
         startvalue++;
         x = y = startvalue;
         i++;
      }
   }
   return d;
}


#endif
