#include <tommath.h>
#ifdef BN_MP_FFT_FLOAT_C

/*
   Even with just 6 bit long chunks it is only good to less
   than ca 80 60-bit limbs (resulting in a 2048 chunks large float array)
   It is also always slower in that range
*/


#include <math.h>
/* 60 bit long limbs with 64 bit machines, so we need a tenth (6 bits) */
#define MP_DIGIT_SIZE (1L<<DIGIT_BIT)
//#define MP_DIGIT_BIT_QUARTER (DIGIT_BIT>>2)
//#define MP_DIGIT_QUARTER (1L<< MP_DIGIT_BIT_QUARTER )
//#define MP_DIGIT_MASK (MP_DIGIT_QUARTER-1)

//#define MP_DIGIT_BIT_FIFTH (DIGIT_BIT/5)
//#define MP_DIGIT_FIFTH (1L<< MP_DIGIT_BIT_FIFTH )
//#define MP_DIGIT_MASK (MP_DIGIT_FIFTH-1)

#define MP_DIGIT_BIT_TENTH (DIGIT_BIT/10)
#define MP_DIGIT_TENTH (1L<< MP_DIGIT_BIT_TENTH )
#define MP_DIGIT_MASK (MP_DIGIT_TENTH-1)


/* base two integer logarithm */
static int highbit(int n)
{
   int r=0;
   int m=n;
   while (m >>= 1) {
      r++;
   }
   return r;
}


/* Transform multiplicands into floating point numbers with TENTH sized digits*/
int mp_dp_to_fft_float(mp_int *a, float **fa,
                 mp_int *b, float **fb, int *length)
{
   int length_a, length_b, length_needed, i, hb, rest;
   float *fft_array_a,*fft_array_b;

   /* Check of the multiplicands happens earlier */
   length_a = a->used;
   length_b = b->used;

   /* Digits get split in TENTHs, so five times the length is needed*/
   length_needed = ((length_a + length_b ))*10 ;
   /* final length must be a power of two to keep the FFTs simple */
   hb = highbit((unsigned  long) length_needed);
   /* check for the rare case that it is already a power of 2 */
   if (length_needed != 1<<hb) {
      length_needed = 1<<(hb+1);
   }
   //fprintf(stderr,"length_needed %d\n",length_needed );
   /* Send computed length back to caller */
   *length = length_needed;

   fft_array_a = XMALLOC(sizeof(float) * (length_needed + 10));
   if (fft_array_a == NULL) {
      return MP_MEM;
   }


   fft_array_b = XMALLOC(sizeof(float) * (length_needed + 10));
   if (fft_array_b == NULL) {
      return MP_MEM;
   }

   for (i = 0; i<length_needed/10; i++) {
      if (i < length_a) {
         fft_array_a[(10*i)]   = (float) (a->dp[i]                            & MP_DIGIT_MASK);
         fft_array_a[(10*i)+1] = (float)((a->dp[i] >>    MP_DIGIT_BIT_TENTH)  & MP_DIGIT_MASK);
         fft_array_a[(10*i)+2] = (float)((a->dp[i] >> (2*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_a[(10*i)+3] = (float)((a->dp[i] >> (3*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_a[(10*i)+4] = (float)((a->dp[i] >> (4*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_a[(10*i)+5] = (float)((a->dp[i] >> (5*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_a[(10*i)+6] = (float)((a->dp[i] >> (6*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_a[(10*i)+7] = (float)((a->dp[i] >> (7*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_a[(10*i)+8] = (float)((a->dp[i] >> (8*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_a[(10*i)+9] = (float)((a->dp[i] >> (9*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
      }
      /* padding a */
      if (i >= length_a) {
         fft_array_a[(10*i)]   = 0.0;
         fft_array_a[(10*i)+1] = 0.0;
         fft_array_a[(10*i)+2] = 0.0;
         fft_array_a[(10*i)+3] = 0.0;
         fft_array_a[(10*i)+4] = 0.0;
         fft_array_a[(10*i)+5] = 0.0;
         fft_array_a[(10*i)+6] = 0.0;
         fft_array_a[(10*i)+7] = 0.0;
         fft_array_a[(10*i)+8] = 0.0;
         fft_array_a[(10*i)+9] = 0.0;

      }
      if (i < length_b) {
         fft_array_b[(10*i)]   = (float) (b->dp[i]                            & MP_DIGIT_MASK);
         fft_array_b[(10*i)+1] = (float)((b->dp[i] >>    MP_DIGIT_BIT_TENTH)  & MP_DIGIT_MASK);
         fft_array_b[(10*i)+2] = (float)((b->dp[i] >> (2*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_b[(10*i)+3] = (float)((b->dp[i] >> (3*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_b[(10*i)+4] = (float)((b->dp[i] >> (4*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_b[(10*i)+5] = (float)((b->dp[i] >> (5*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_b[(10*i)+6] = (float)((b->dp[i] >> (6*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_b[(10*i)+7] = (float)((b->dp[i] >> (7*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_b[(10*i)+8] = (float)((b->dp[i] >> (8*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_b[(10*i)+9] = (float)((b->dp[i] >> (9*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);

      }
      /* padding b */
      if (i >= length_b) {
         fft_array_b[(10*i)]   = 0.0;
         fft_array_b[(10*i)+1] = 0.0;
         fft_array_b[(10*i)+2] = 0.0;
         fft_array_b[(10*i)+3] = 0.0;
         fft_array_b[(10*i)+4] = 0.0;
         fft_array_b[(10*i)+5] = 0.0;
         fft_array_b[(10*i)+6] = 0.0;
         fft_array_b[(10*i)+7] = 0.0;
         fft_array_b[(10*i)+8] = 0.0;
         fft_array_b[(10*i)+9] = 0.0;

      }
   }
   // there is a small problem with divisibility of 2^n and 10, so ...

   rest = (length_needed/10)*10;
   for(i=rest;i<length_needed + 10;i++){
      fft_array_a[i] = 0.0;
      fft_array_b[i] = 0.0;
   }

   /* Send the route to memory back to caller */
   *fa = fft_array_a;
   *fb = fft_array_b;
   return MP_OKAY;
}


/* same as dp_to_fft() for a single multiplicand for squaring */
int mp_dp_to_fft_single_float(mp_int *a, float **fa, int *length)
{
   int length_a,  length_needed, i, hb, rest;
   float *fft_array_a;
   length_a = a->used;
   length_needed = (length_a * 2)*5 ;
   hb = highbit((unsigned  long) length_needed);
   if (length_needed != 1<<hb) {
      length_needed = 1<<(hb+1);
   }
   *length = length_needed;
   fft_array_a = XMALLOC(sizeof(float) * (length_needed + 5));
   if (fft_array_a == NULL) {
      return MP_MEM;
   }
   for (i = 0; i<length_needed/5; i++) {
      if (i < length_a) {
         fft_array_a[(5*i)]   = (float)( a->dp[i]                            & MP_DIGIT_MASK);
         fft_array_a[(5*i)+1] = (float)((a->dp[i] >>    MP_DIGIT_BIT_TENTH)  & MP_DIGIT_MASK);
         fft_array_a[(5*i)+2] = (float)((a->dp[i] >> (2*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_a[(5*i)+3] = (float)((a->dp[i] >> (3*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
         fft_array_a[(5*i)+4] = (float)((a->dp[i] >> (4*MP_DIGIT_BIT_TENTH)) & MP_DIGIT_MASK);
      }
      if (i >= length_a) {
         fft_array_a[(5*i)]   = 0.0;
         fft_array_a[(5*i)+1] = 0.0;
         fft_array_a[(5*i)+2] = 0.0;
         fft_array_a[(5*i)+3] = 0.0;
         fft_array_a[(5*i)+4] = 0.0;
      }
   }
   rest = (length_needed/5)*5;
   for(i=rest;i<length_needed + 5;i++){
      fft_array_a[i] = 0.0;
   }
   *fa = fft_array_a;
   return MP_OKAY;
}

int mp_fft_to_dp_float(float *fft_array, mp_int *a,int length)
{
   int new_length, i,j,e;
   mp_word carry = 0,temp;

   /* Result cannot exceed length/2, hence add two */
   new_length = length;

   /* Preallocate some memory for the result. */
   if (a->alloc < new_length) {
      if ((e = mp_grow(a, new_length)) != MP_OKAY) {
         return e;
      }
   }

   /* The FFT multiplication does no carry (it's one of the tricks of it) */

   /* Hard to paralellize because of the carry */
   for (i=0; i<length; i++) {
      temp = carry;
      carry = 0;
      temp  += (mp_word)(roundf(fft_array[i]));
      if (temp >= MP_DIGIT_TENTH) {
         carry = temp / (mp_word)MP_DIGIT_TENTH;
         temp  = temp % (mp_word)MP_DIGIT_TENTH;
      }
      /* memory is still expensive, not a thing to waste easily */
      fft_array[i] = (float)temp;
   }

#if __STDC_VERSION__ >= 199901L
#define NEEDS_FE_RESET 1
#include <fenv.h>
   fenv_t envp;
   /* backup of floating point environment settings */
   if (fegetenv(&envp)) return MP_VAL;
   /* Set rounding mode to "nearest". Default, but better safe than sorry */
   if (fesetround(FE_TONEAREST)) return MP_VAL;
#endif


   /* re-marry the digits */
   for (i=0,j=0; j<new_length; i++,j+=10) {
      a->dp[i]  = (mp_digit)(roundf(fft_array[j+9]))  & MP_DIGIT_MASK;
      a->dp[i] <<= MP_DIGIT_BIT_TENTH;
      a->dp[i] |= (mp_digit)(roundf(fft_array[j+8]))  & MP_DIGIT_MASK;
      a->dp[i] <<= MP_DIGIT_BIT_TENTH;
      a->dp[i] |= (mp_digit)(roundf(fft_array[j+7]))  & MP_DIGIT_MASK;
      a->dp[i] <<= MP_DIGIT_BIT_TENTH;
      a->dp[i] |= (mp_digit)(roundf(fft_array[j+6]))  & MP_DIGIT_MASK;
      a->dp[i] <<= MP_DIGIT_BIT_TENTH;
      a->dp[i] |= (mp_digit)(roundf(fft_array[j+5]))  & MP_DIGIT_MASK;
      a->dp[i] <<= MP_DIGIT_BIT_TENTH;
      a->dp[i] |= (mp_digit)(roundf(fft_array[j+4]))  & MP_DIGIT_MASK;
      a->dp[i] <<= MP_DIGIT_BIT_TENTH;
      a->dp[i] |= (mp_digit)(roundf(fft_array[j+3]))  & MP_DIGIT_MASK;
      a->dp[i] <<= MP_DIGIT_BIT_TENTH;
      a->dp[i] |= (mp_digit)(roundf(fft_array[j+2]))  & MP_DIGIT_MASK;
      a->dp[i] <<= MP_DIGIT_BIT_TENTH;
      a->dp[i] |= (mp_digit)(roundf(fft_array[j+1]))  & MP_DIGIT_MASK;
      a->dp[i] <<= MP_DIGIT_BIT_TENTH;
      a->dp[i] |= (mp_digit)(roundf(fft_array[j]))    & MP_DIGIT_MASK;
      /* and count them all */
      a->used++;
   }
  // fprintf(stderr,"a.used %d\n",a->used );
   if (carry) {
      a->dp[i] = carry;
      a->used++;
   }
   mp_clamp(a);
//   fprintf(stderr,"a.used %d\n",a->used );
   return MP_OKAY;
}



/*
  The size of the L1-cache in bytes. The number here is that of the data cache
  part of an AMD Duron. The Linux kernel gives a lot of information e.g.:
    grep . /sys/devices/system/cpu/cpu0/cache/index*//*
  There is also lscpu(1) wich is easier to use.
  On Windows:
    http://msdn.microsoft.com/en-us/library/ms683194.aspx
    http://www.cpuid.com/softwares/cpu-z.htm
  Lack of access to a Mac leaves that part blank. The new MacOS is based on BSD,
  so 'dmesg' might work or
    cat /var/run/dmesg.boot | grep CPU
 */
#ifndef L1_SIZE
//#define L1_SIZE 65536
#define L1_SIZE 16384
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937511
#endif
#define TWOPI (2.0*M_PI)


static void fht_dif_iterative_float(float *x, unsigned long n, int do_loop)
{
   unsigned long m,mh,mq;
   unsigned long i,j,k;
   float a,b,t, c,s, u,v,tmp;
   float *dp;
   for (m=n; m > 1; m >>= 1) {
      mh = m >> 1;
      mq = mh >> 1;
      t = M_PI / (float)mh;
      a = sinf(0.5 * t);
      a *= 2.0 * a;
      b = sinf(t);
      for (i = 0; i < n; i += m) {
         dp = x + i;
         for (j = 0, k = mh; j < mh; ++j, ++k) {
            u = dp[j];
            v = dp[k];
            dp[j] = u + v;
            dp[k] = u - v;
         }
         dp += mh;
         c = 1.0;
         s = 0.0;
         for (j = 1, k = mh - 1; j < mq; ++j, --k) {
            tmp = c;
            c -= a * c + b * s;
            s -= a * s - b * tmp;
            u = dp[j];
            v = dp[k];
            dp[j] = u * c + v * s;
            dp[k] = u * s - v * c;
         }
      }
      if (!do_loop)break;
   }
   return;
}



static void fht_dif_rec_float(float *x, unsigned long n)
{
   unsigned long nh;
   if (n == 1)
      return;
   if (n < (unsigned long)(L1_SIZE / (2 * sizeof(float)))) {
      fht_dif_iterative_float(x, n, 1);
      return;
   }
   fht_dif_iterative_float(x, n, 0);
   nh = n >> 1;
   fht_dif_rec_float(x, nh);
   fht_dif_rec_float(x + nh, nh);
   return;
}

static void fht_dit_iterative_float(float *x, unsigned long n, int do_loop)
{
   unsigned long m, mh ,mq;
   unsigned long i,j,k;
   float a,b,t, u,v, c,s, tmp;
   float *dp;

   m = (do_loop)?2:n;
   for (; m <= n; m <<= 1) {
      mh = m >> 1;
      mq = mh >> 1;
      t = M_PI / (float)mh;
      a = sinf(0.5 * t);
      a *= 2.0 * a;
      b = sinf(t);
      for (i = 0; i < n; i += m) {
         dp = x + i + mh;
         c = 1.0;
         s = 0.0;
         for (j = 1, k = mh - 1; j < mq; ++j, --k) {
            tmp = c;
            c -= a * c + b * s;
            s -= a * s - b * tmp;
            u = dp[j];
            v = dp[k];
            dp[j] = u * c + v * s;
            dp[k] = u * s - v * c;
         }
         dp -= mh;
         for (j = 0, k = mh; j < mh; ++j, ++k) {
            u = dp[j];
            v = dp[k];
            dp[j] = u + v;
            dp[k] = u - v;
         }
      }
   }
   return;
}

static void fht_dit_rec_float(float *x, unsigned long n)
{
   unsigned long nh;

   if (n == 1)
      return;
   if (n < (unsigned long)(L1_SIZE / (2 * sizeof(float)))) {
      fht_dit_iterative_float(x,n,1);
      return;
   }
   nh = n >> 1;
   fht_dit_rec_float(x, nh);
   fht_dit_rec_float(x + nh, nh);
   fht_dit_iterative_float(x,n,0);
   return;
   return;
}


static void fht_conv_core_float(float *f, float *g,unsigned long n, float v/*=0.0*/)
{
   unsigned long nh,r,rm,k,km,tr,m;
   float xi,xj, yi,yj;
   if (v==0.0)  v = 1.0/n;

   g[0] *= (v * f[0]);
   if (n>=2) g[1] *= (v * f[1]);
   if (n<4)   return;
   v *= 0.5;
   nh = (n>>1);
   r=nh;
   rm=n-1;
   xi = f[r];
   xj = f[rm];
   yi = g[r];
   yj = g[rm];
   g[r]  = v*((xi + xj)*yi + (xi - xj)*yj);
   g[rm] = v*((-xi + xj)*yi + (xi + xj)*yj);

   k=2;
   km=n-2;
   while (k<nh) {
      rm -= nh;
      tr = r;
      r^=nh;
      for (m=(nh>>1); !((r^=m)&m); m>>=1) {
         ;
      }

      xi = f[r];
      xj = f[rm];
      yi = g[r];
      yj = g[rm];
      g[r]  = v*((xi + xj)*yi + (xi - xj)*yj);
      g[rm] = v*((-xi + xj)*yi + (xi + xj)*yj);
      --km;
      ++k;
      rm += (tr-r);
      r += nh;
      xi = f[r];
      xj = f[rm];
      yi = g[r];
      yj = g[rm];
      g[r]  = v*((xi + xj)*yi + (xi - xj)*yj);
      g[rm] = v*((-xi + xj)*yi + (xi + xj)*yj);
      --km;
      ++k;
   }
   return;

}

static void fht_autoconv_core_float(float *f,unsigned long n, float v/*=0.0*/)
{
   unsigned long nh,r,rm,k,km,tr,m;
   float xi,xj, xi2, xj2,xij ;
   if (v==0.0)  v = 1.0/n;

   f[0] *= (v * f[0]);
   if (n>=2) f[1] *= (v * f[1]);
   if (n<4)   return;
   v *= 0.5;
   nh = (n>>1);
   r=nh;
   rm=n-1;
   xi = f[r];
   xj = f[rm];
   xi2 = xi*xi;
   xj2 = xj*xj;
   xij = (2*xi*xj);
   f[r]  = v*(xi2 + xij - xj2);
   f[rm] = v*(-xi2 + xij + xj2);

   k=2;
   km=n-2;
   while (k<nh) {
      rm -= nh;
      tr = r;
      r^=nh;
      for (m=(nh>>1); !((r^=m)&m); m>>=1) {
         ;
      }
      xi = f[r];
      xj = f[rm];
      xi2 = xi*xi;
      xj2 = xj*xj;
      xij = (2*xi*xj);
      f[r]  = v*(xi2 + xij - xj2);
      f[rm] = v*(-xi2 + xij + xj2);

      --km;
      ++k;
      rm += (tr-r);
      r += nh;
      xi = f[r];
      xj = f[rm];
      xi2 = xi*xi;
      xj2 = xj*xj;
      xij = (2*xi*xj);
      f[r]  = v*(xi2 + xij - xj2);
      f[rm] = v*(-xi2 + xij + xj2);
      --km;
      ++k;
   }
   return;
}



/* Public: FHT convolution */
int mp_fft_float(float *x, float *y, unsigned long length)
{
   unsigned long n;
   n = (length);
   if (n < 2) return MP_VAL;
   fht_dif_rec_float(x,(n));
   fht_dif_rec_float(y,(n));
   fht_conv_core_float(x, y,(n), 0.0);
   fht_dit_rec_float(y, (n));
   return MP_OKAY;
}
/* Public: FHT auto-convolution */
int mp_fft_sqr_d_float(float *x, unsigned long length)
{
   unsigned long n;
   n = (length);
   if (n < 2) return MP_VAL;
   fht_dif_rec_float(x,(n));
   fht_autoconv_core_float(x,(n), 0.0);
   fht_dit_rec_float(x, (n));
   return MP_OKAY;
}


#if (__STDC_VERSION__ >= 199901L) &&  (NEEDS_FE_RESET == 1)
/* Reset floating point environment settings  */
if (fesetenv(envp)) return MP_VAL;
#endif


#endif
