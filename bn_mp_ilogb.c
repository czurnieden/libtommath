#include <tommath.h>
#ifdef BN_MP_ILOGB_C

int mp_ilogb(mp_int * a, mp_int * base, mp_int * c) {
  int err, bits, cmp;
  mp_int low, bracket_low, high, bracket_high, mid, bracket_mid, t;

  err = MP_OKAY;
  if (a->sign == MP_NEG || base->sign == MP_NEG) {
    return MP_VAL;
  }
  if (base->used <= 1) {
    if (base->dp[0] == 2) {
      bits = mp_count_bits(a) - 1;
      mp_set(c, bits);
      return err;
    }
    if (base->dp[0] < 2) {
      return MP_VAL;
    }
  }
  if (mp_cmp_d(a, (mp_digit) (0)) == MP_EQ) {
    return MP_VAL;
  }
  if (mp_cmp_d(a, (mp_digit) (1)) == MP_EQ) {
    mp_zero(c);
    return err;
  }
  cmp = mp_cmp(a, base);
  if (cmp == MP_LT) {
    mp_zero(c);
    return err;
  } else if (cmp == MP_EQ) {
    mp_set(c, (mp_digit) (1));
    return err;
  }

  if ((err =
       mp_init_multi(&low, &bracket_low, &high, &bracket_high, &mid,
		     &bracket_mid, &t, NULL)) != MP_OKAY) {
    return err;
  }
  mp_zero(&low);
  mp_set(&bracket_low, 1);
  mp_set(&high, 1);
  if ((err = mp_copy(base, &bracket_high)) != MP_OKAY) {
    goto _ERR;
  }

  while (mp_cmp(&bracket_high, a) == MP_LT) {
    if ((err = mp_copy(&high, &low)) != MP_OKAY) {
      goto _ERR;
    }
    if ((err = mp_copy(&bracket_high, &bracket_low)) != MP_OKAY) {
      goto _ERR;
    }
    if ((err = mp_mul_2d(&high, 1, &high)) != MP_OKAY) {
      goto _ERR;
    }
    if ((err = mp_sqr(&bracket_high, &bracket_high)) != MP_OKAY) {
      goto _ERR;
    }
  }
  if (mp_cmp(a, &bracket_high) == MP_EQ) {
    mp_exch(&high, c);
    goto _ERR;
  }
  mp_sub(&high, &low, &t);
  while (mp_cmp_d(&t, (mp_digit) (1)) == MP_GT) {
    if ((err = mp_add(&low, &high, &t)) != MP_OKAY) {
      goto _ERR;
    }
    if ((err = mp_div_2d(&t, 1, &mid, NULL)) != MP_OKAY) {
      goto _ERR;
    }
    if ((err = mp_sub(&mid, &low, &t)) != MP_OKAY) {
      goto _ERR;
    }
    if ((err = mp_expt(base, &t, &t)) != MP_OKAY) {
      goto _ERR;
    }
    if ((err = mp_mul(&bracket_low, &t, &bracket_mid)) != MP_OKAY) {
      goto _ERR;
    }
    cmp = mp_cmp(a, &bracket_mid);
    if (cmp == MP_LT) {
      mp_exch(&mid, &high);
      mp_exch(&bracket_mid, &bracket_high);
    }
    if (cmp == MP_GT) {
      mp_exch(&mid, &low);
      mp_exch(&bracket_mid, &bracket_low);
    }
    if (cmp == MP_EQ) {
      mp_exch(&mid, c);
      goto _ERR;
    }
    if ((err = mp_sub(&high, &low, &t)) != MP_OKAY) {
      goto _ERR;
    }
  }

  if (mp_cmp(&bracket_high, a) == MP_EQ) {
    mp_exch(&high, c);
    goto _ERR;
  } else {
    mp_exch(&low, c);
    goto _ERR;
  }
_ERR:
  mp_clear_multi(&low, &bracket_low, &high, &bracket_high, &mid, &bracket_mid,
		 &t, NULL);
  return err;
}

#endif

