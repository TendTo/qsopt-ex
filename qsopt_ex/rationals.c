/*
 * This file is part of the Yices SMT Solver.
 * Copyright (C) 2017 SRI International.
 *
 * Yices is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Yices is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Yices.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Operations on rational numbers
 * - rationals are represented as pairs of 32 bit integers
 *   or if they are too large as a pointer to a gmp rational.
 */

#include "rationals.h"

#include <assert.h>
#include <gmp.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "eg_mem.h"
#include "logging-private.h"
#include "mpq_aux.h"
#include "mpq_stores.h"

/*
 * gcd of two 32bit unsigned positive numbers.
 */
static inline uint32_t gcd32(uint32_t a, uint32_t b) {
  assert(a > 0 && b > 0);

  uint32_t k = 1;
  uint32_t x = a | b;
  while ((x & 1) == 0) {
    k <<= 1;
    x >>= 1;
    a >>= 1;
    b >>= 1;
  }

  do {
    if ((a & 1) == 0) {
      a >>= 1;
    } else if ((b & 1) == 0) {
      b >>= 1;
    } else if (a >= b) {
      a = (a - b) >> 1;
    } else {
      b = (b - a) >> 1;
    }
  } while (a > 0);

  return b * k;
}

/*
 * gcd of two 64bit unsigned positive numbers
 */
static inline uint64_t gcd64(uint64_t a, uint64_t b) {
  assert(a > 0 && b > 0);

  uint64_t k = 1;
  uint64_t x = a | b;
  while ((x & 1) == 0) {
    k <<= 1;
    x >>= 1;
    a >>= 1;
    b >>= 1;
  }

  do {
    if ((a & 1) == 0) {
      a >>= 1;
    } else if ((b & 1) == 0) {
      b >>= 1;
    } else if (a >= b) {
      a = (a - b) >> 1;
    } else {
      b = (b - a) >> 1;
    }
  } while (a > 0);

  return b * k;
}

static mpq_store_t mpq_store;

/*
 *  String buffer for parsing.
 */
static char* string_buffer = NULL;
static uint32_t string_buffer_length = 0;

/*
 * Print an error then abort on division by zero
 */
static void division_by_zero(void) {
  fprintf(stderr, "\nRationals: division by zero\n");
  abort();
}

/*
 * Initialize everything including the string lock
 * if we're in thread-safe mode.
 */
void init_rationals(void) {
  init_mpq_aux();
  init_mpqstore(&mpq_store);
  string_buffer = NULL;
  string_buffer_length = 0;
}

/*
 * Cleanup: free memory
 */
void cleanup_rationals(void) {
  cleanup_mpq_aux();
  delete_mpqstore(&mpq_store);
  EGmalloc(string_buffer);
}

/*************************
 *  MPQ ALLOCATION/FREE  *
 ************************/

/*
 * Allocates a new mpq object.
 */
static inline mpq_ptr new_mpq(void) { return mpqstore_alloc(&mpq_store); }

/*
 * Deallocates a new mpq object
 */
static void release_mpq(rat_t r) {
  assert(rat_is_gmp(r));
  mpqstore_free(&mpq_store, rat_get_mpq(r));
}

/*
 * Free mpq number attached to r if any, then set r to 0/1.
 * Must be called before r is deleted to prevent memory leaks.
 */
void rat_clear(rat_t r) {
  if (rat_is_gmp(r)) {
    release_mpq(r);
  }
  r->s.num = 0;
  r->s.den = RAT_ONE_DEN;
}

/*******************
 *  NORMALIZATION  *
 ******************/

/*
 * Bounds on numerator and denominators.
 *
 * a/b can be safely stored as a pair (int32_t/uint32_t)
 * if MIN_NUMERATOR <= a <= MAX_NUMERATOR
 * and 1 <= b <= MAX_DENOMINATOR
 *
 * otherwise a/b must be stored as a gmp rational.
 *
 * The bounds are such that
 * - (a/1)+(b/1), (a/1) - (b/1) can be computed using
 *   32bit arithmetic without overflow.
 * - a/b stored as a pair implies -a/b and b/a can be stored
 *   as pairs too.
 */
#define MAX_NUMERATOR (INT32_MAX >> 1)
#define MIN_NUMERATOR (-MAX_NUMERATOR)
#define MAX_DENOMINATOR MAX_NUMERATOR

/*
 * Store num/den in r
 */
static inline void set_rat32(rat_t r, int32_t num, uint32_t den) {
  assert(MIN_NUMERATOR <= num && num <= MAX_NUMERATOR && den <= MAX_DENOMINATOR);
  r->s.den = den << 1;
  r->s.num = num;
}

/*
 * Normalization: construct rational a/b when
 * a and b are two 64bit numbers.
 * - b must be non-zero
 */
void rat_set_int64(rat_t r, int64_t a, uint64_t b) {
  uint64_t abs_a;
  mpq_ptr q;
  bool a_positive;

  assert(b > 0);

  if (a == 0 || (b == 1 && MIN_NUMERATOR <= a && a <= MAX_NUMERATOR)) {
    if (rat_is_gmp(r)) {
      release_mpq(r);
    }
    set_rat32(r, a, 1);
    return;
  }

  // absolute value and sign of a.
  if (a >= 0) {
    abs_a = (uint64_t)a;
    a_positive = true;
  } else {
    abs_a = (uint64_t)-a;  // Note: this works even when a = -2^63
    a_positive = false;
  }

  // abs_a and b are positive. remove powers of 2
  while (((abs_a | b) & 15) == 0) {
    abs_a >>= 4;
    b >>= 4;
  }
  switch ((abs_a | b) & 7) {
    case 0:
      abs_a >>= 3;
      b >>= 3;
      break;
    case 1:
      break;
    case 2:
      abs_a >>= 1;
      b >>= 1;
      break;
    case 3:
      break;
    case 4:
      abs_a >>= 2;
      b >>= 2;
      break;
    case 5:
      break;
    case 6:
      abs_a >>= 1;
      b >>= 1;
      break;
    case 7:
      break;
  }

  // abs_a and b are positive, and at least one of them is odd.
  // if abs_a <= 2 or b <= 2 then gcd = 1.
  if (abs_a > 2 && b > 2) {
    uint64_t a_1 = abs_a;
    uint64_t b_1 = b;

    // compute gcd of abs_a and b
    // loop invariant: abs_a is odd or b is odd (or both)
    for (;;) {
      if ((a_1 & 1) == 0) {
        a_1 >>= 1;
      } else if ((b_1 & 1) == 0) {
        b_1 >>= 1;
      } else if (a_1 >= b_1) {
        a_1 = (a_1 - b_1) >> 1;
        if (a_1 == 0) break;
      } else {
        b_1 = (b_1 - a_1) >> 1;
      }
    }

    // b_1 is gcd(abs_a, b)
    if (b_1 != 1) {
      abs_a /= b_1;
      b /= b_1;
    }
  }

  // abs_a and b are mutually prime and positive

  // restore a
  a = a_positive ? ((int64_t)abs_a) : -((int64_t)abs_a);

  // assign to r
  if (abs_a <= MAX_NUMERATOR && b <= MAX_DENOMINATOR) {
    if (rat_is_gmp(r)) {
      release_mpq(r);
    }
    set_rat32(r, (int32_t)a, (uint32_t)b);
  } else {
    if (rat_is_gmp(r)) {
      q = rat_get_mpq(r);
    } else {
      q = new_mpq();
      rat_set_gmp(r, q);
    }
    mpq_set_int64(q, a, b);
  }
}

/*
 * Normalization: construct a/b when a and b are 32 bits
 */
void rat_set_int32(rat_t r, int32_t a, uint32_t b) {
  uint32_t abs_a;
  mpq_ptr q;
  bool a_positive;

  assert(b > 0);

  if (a == 0 || (b == 1 && MIN_NUMERATOR <= a && a <= MAX_NUMERATOR)) {
    if (rat_is_gmp(r)) {
      release_mpq(r);
    }
    set_rat32(r, a, 1);
    return;
  }

  // absolute value and sign of a.
  if (a >= 0) {
    abs_a = (uint32_t)a;
    a_positive = true;
  } else {
    abs_a = (uint32_t)-a;  // Note: this works even when a = -2^31
    a_positive = false;
  }

  // abs_a and b are positive. remove powers of 2
  while (((abs_a | b) & 15) == 0) {
    abs_a >>= 4;
    b >>= 4;
  }
  switch ((abs_a | b) & 7) {
    case 0:
      abs_a >>= 3;
      b >>= 3;
      break;
    case 1:
      break;
    case 2:
      abs_a >>= 1;
      b >>= 1;
      break;
    case 3:
      break;
    case 4:
      abs_a >>= 2;
      b >>= 2;
      break;
    case 5:
      break;
    case 6:
      abs_a >>= 1;
      b >>= 1;
      break;
    case 7:
      break;
  }

  // abs_a and b are positive, and at least one of them is odd.
  // if abs_a <= 2 or b <= 2 then gcd = 1.
  if (abs_a > 2 && b > 2) {
    uint32_t a_1 = abs_a;
    uint32_t b_1 = b;

    // compute gcd of abs_a and b
    // loop invariant: abs_a is odd or b is odd (or both)
    for (;;) {
      if ((a_1 & 1) == 0) {
        a_1 >>= 1;
      } else if ((b_1 & 1) == 0) {
        b_1 >>= 1;
      } else if (a_1 >= b_1) {
        a_1 = (a_1 - b_1) >> 1;
        if (a_1 == 0) break;

      } else {
        b_1 = (b_1 - a_1) >> 1;
      }
    }

    // b_1 is gcd(abs_a, b)
    if (b_1 != 1) {
      abs_a /= b_1;
      b /= b_1;
    }
  }

  // abs_a and b are mutually prime and positive

  // restore a
  a = a_positive ? ((int32_t)abs_a) : -((int32_t)abs_a);

  // assign to r
  if (abs_a <= MAX_NUMERATOR && b <= MAX_DENOMINATOR) {
    if (rat_is_gmp(r)) {
      release_mpq(r);
    }
    set_rat32(r, (int32_t)a, (uint32_t)b);
  } else {
    if (rat_is_gmp(r)) {
      q = rat_get_mpq(r);
    } else {
      q = new_mpq();
      rat_set_gmp(r, q);
    }
    mpq_set_int32(q, a, b);
  }
}

/*
 * Construct r = a/1
 */
void rat_set64(rat_t r, int64_t a) {
  mpq_ptr q;

  if (MIN_NUMERATOR <= a && a <= MAX_NUMERATOR) {
    if (rat_is_gmp(r)) {
      release_mpq(r);
    }
    set_rat32(r, (int32_t)a, 1);
  } else {
    if (!rat_is_gmp(r)) {
      q = new_mpq();
      rat_set_gmp(r, q);
    } else {
      q = rat_get_mpq(r);
    }
    mpq_set_int64(q, a, 1);
  }
}

void rat_set32(rat_t r, int32_t a) {
  mpq_ptr q;

  if (MIN_NUMERATOR <= a && a <= MAX_NUMERATOR) {
    if (rat_is_gmp(r)) {
      release_mpq(r);
    }
    set_rat32(r, a, 1);
  } else {
    if (!rat_is_gmp(r)) {
      q = new_mpq();
      rat_set_gmp(r, q);
    } else {
      q = rat_get_mpq(r);
    }
    mpq_set_int32(q, a, 1);
  }
}

/*
 * Convert r to a gmp number.
 * r->num and r->den must have no common factor
 * and r->den must be non-zero.
 */
static void convert_to_gmp(rat_t r) {
  mpq_ptr q;

  assert(!rat_is_gmp(r));

  q = new_mpq();
  mpq_set_int32(q, rat_get_num(r), rat_get_den(r));
  rat_set_gmp(r, q);
}

/*
 * Set r to a gmp with value = a.
 */
static void set_to_gmp64(rat_t r, int64_t a) {
  mpq_ptr q;

  assert(!rat_is_gmp(r));

  q = new_mpq();
  mpq_set_int64(q, a, 1);
  rat_set_gmp(r, q);
}

/*****************
 *  ASSIGNMENTS  *
 ****************/

/*
 * Convert mpq to a pair of integers if possible.
 */
void rat_normalize(rat_t r) {
  mpq_ptr q;
  unsigned long den;
  long num;

  if (rat_is_gmp(r)) {
    q = rat_get_mpq(r);
    if (mpz_fits_ulong_p(mpq_denref(q)) && mpz_fits_slong_p(mpq_numref(q))) {
      num = mpz_get_si(mpq_numref(q));
      den = mpz_get_ui(mpq_denref(q));
      if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR && den <= MAX_DENOMINATOR) {
        mpqstore_free(&mpq_store, q);
        set_rat32(r, (int32_t)num, (uint32_t)den);
      }
    }
  }
}

static void rat_denormalize(rat_t r) {
  mpq_ptr q;

  if (rat_is_rat32(r)) {
    q = new_mpq();
    mpq_set_si(q, rat_get_num(r), rat_get_den(r));
    rat_set_gmp(r, q);
  }
  assert(rat_is_gmp(r));
}

/*
 * Prepare to assign an mpq number to r
 * - if r is not a mpq number, allocate it
 */
static inline void rat_prepare(rat_t r) {
  if (!rat_is_gmp(r)) {
    rat_set_gmp(r, new_mpq());
  }
}

/*
 * assign r:= z/1
 */
void rat_set_mpz(rat_t r, const mpz_t z) {
  mpq_ptr q;

  rat_prepare(r);
  q = rat_get_mpq(r);
  mpq_set_z(q, z);
  rat_normalize(r);
}

void rat_set_mpf(rat_t r, const mpf_t f) {
  mpq_ptr q;

  rat_prepare(r);
  q = rat_get_mpq(r);
  mpq_set_f(q, f);
  rat_normalize(r);
}

/*
 * Copy q into r
 */
void rat_set_mpq(rat_t r, const mpq_t q) {
  mpq_ptr qt;

  rat_prepare(r);
  qt = rat_get_mpq(r);
  mpq_set(qt, q);
  rat_normalize(r);
}

void rat_set_mpq_prenormalize(rat_t r, const mpq_t q) {
  unsigned long den;
  long num;
  mpq_ptr qt;

  // If the numerator and denominator of q fit in 32 bits,
  // we can store it as a pair of integers
  if (mpz_fits_ulong_p(mpq_denref(q)) && mpz_fits_slong_p(mpq_numref(q))) {
    num = mpz_get_si(mpq_numref(q));
    den = mpz_get_ui(mpq_denref(q));
    if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR && den <= MAX_DENOMINATOR) {
      set_rat32(r, (int32_t)num, (uint32_t)den);
      assert(rat_is_rat32(r));
      return;
    }
  }

  // Otherwise, we need to allocate a new mpq number and copy q into it
  rat_prepare(r);
  qt = rat_get_mpq(r);
  mpq_set(qt, q);
  assert(rat_is_gmp(r));
}

/*
 * Sets a rational to a value from a double
 */
void rat_set_double(rat_t r, double d) {
  mpq_ptr qt;

  rat_prepare(r);
  qt = rat_get_mpq(r);
  mpq_set_d(qt, d);
  rat_normalize(r);
}

/*
 * Copy r2 into r1
 */
void rat_set(rat_t r1, const rat_t r2) {
  mpq_ptr q1, q2;

  if (rat_is_gmp(r2)) {
    rat_prepare(r1);
    q1 = rat_get_mpq(r1);
    q2 = rat_get_mpq(r2);
    mpq_set(q1, q2);
  } else {
    if (rat_is_gmp(r1)) {
      release_mpq(r1);
    }
    r1->s.num = r2->s.num;
    r1->s.den = r2->s.den;
  }
}

/*
 * Copy opposite of r2 into r1
 */
void rat_set_neg(rat_t r1, const rat_t r2) {
  mpq_ptr q1, q2;

  if (rat_is_gmp(r2)) {
    rat_prepare(r1);
    q1 = rat_get_mpq(r1);
    q2 = rat_get_mpq(r2);
    mpq_neg(q1, q2);
  } else {
    if (rat_is_gmp(r1)) {
      release_mpq(r1);
    }
    r1->s.num = -r2->s.num;
    r1->s.den = r2->s.den;
  }
}

/*
 * Copy the absolute value of r2 into r1
 */
void rat_set_abs(rat_t r1, const rat_t r2) {
  if (rat_is_gmp(r2)) {
    mpq_ptr q1;
    mpq_ptr q2;
    rat_prepare(r1);
    q1 = rat_get_mpq(r1);
    q2 = rat_get_mpq(r2);
    mpq_abs(q1, q2);
  } else {
    if (rat_is_gmp(r1)) {
      release_mpq(r1);
    }
    r1->s.den = r2->s.den;
    if (r2->s.num < 0) {
      r1->s.num = -r2->s.num;
    } else {
      r1->s.num = r2->s.num;
    }
  }
}

/*
 * Copy the numerator of r2 into r1
 * - r1 must be initialized
 * - r2 and r1 must be different objects
 */
void rat_set_num(rat_t r1, const rat_t r2) {
  mpq_ptr q1, q2;
  long num;

  if (rat_is_gmp(r2)) {
    q2 = rat_get_mpq(r2);
    if (mpz_fits_slong_p(mpq_numref(q2))) {
      num = mpz_get_si(mpq_numref(q2));
      if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR) {
        if (rat_is_gmp(r1)) {
          release_mpq(r1);
        }
        set_rat32(r1, num, 1);
        return;
      }
    }
    rat_prepare(r1);
    q1 = rat_get_mpq(r1);
    mpq_set_z(q1, mpq_numref(q2));
  } else {
    if (rat_is_gmp(r1)) {
      release_mpq(r1);
    }
    set_rat32(r1, rat_get_num(r2), 1);
  }
}

// /*
//  * Copy the denominator of r2 into r1
//  * - r1 must be initialized
//  * - r1 and r2 must be different objects
//  */
void rat_set_den(rat_t r1, const rat_t r2) {
  mpq_ptr q1, q2;
  unsigned long den;

  if (rat_is_gmp(r2)) {
    q2 = rat_get_mpq(r2);
    if (mpz_fits_ulong_p(mpq_denref(q2))) {
      den = mpz_get_ui(mpq_denref(q2));
      if (den <= MAX_DENOMINATOR) {
        if (rat_is_gmp(r1)) {
          release_mpq(r1);
        }
        r1->s.num = den;
        r1->s.den = RAT_ONE_DEN;
        return;
      }
    }
    rat_prepare(r1);
    q1 = rat_get_mpq(r1);
    mpq_set_z(q1, mpq_denref(q2));

  } else {
    if (rat_is_gmp(r1)) {
      release_mpq(r1);
    }
    r1->s.num = rat_get_den(r2);
    r1->s.den = RAT_ONE_DEN;
  }
}

/*******************************************
 *  PARSING: CONVERT STRINGS TO RATIONALS  *
 ******************************************/

/*
 * Resize string_buffer if necessary
 */
static void resize_string_buffer(uint32_t new_size) {
  uint32_t n;

  if (string_buffer_length < new_size) {
    /*
     * try to make buffer 50% larger
     * in principle the n += n >> 1 could overflow but
     * then (n < new_size) will be true
     * so n will be fixed to new_size.
     *
     * all this is very unlikely to happen in practice
     * (this would require extremely large strings).
     */
    n = string_buffer_length + 1;
    n += n >> 1;
    if (n < new_size) n = new_size;

    string_buffer = (char*)EGrealloc(string_buffer, n);
    string_buffer_length = n;
  }
}

/*
 * Assign q0 to r and try to convert to a pair of integers.
 * - q0 must be canonicalized
 * - return 0.
 */
static int rat_set_q0(rat_t r, mpq_t* q0) {
  mpq_ptr q;
  unsigned long den;
  long num;

  // try to store q0 as a pair num/den
  if (mpz_fits_ulong_p(mpq_denref(*q0)) && mpz_fits_slong_p(mpq_numref(*q0))) {
    num = mpz_get_si(mpq_numref(*q0));
    den = mpz_get_ui(mpq_denref(*q0));
    if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR && den <= MAX_DENOMINATOR) {
      if (rat_is_gmp(r)) {
        release_mpq(r);
      }
      set_rat32(r, (int32_t)num, (uint32_t)den);
      return 0;
    }
  }

  // copy q0
  if (!rat_is_gmp(r)) {
    rat_prepare(r);
  }
  q = rat_get_mpq(r);
  mpq_set(q, *q0);
  return 0;
}

/*
 * Conversion from a string in format supported by gmp:
 *   <+/->numerator/denominator or <+/->numerator
 * with numerator and denominator in base 10.
 * - returns -1 and leaves r unchanged if s is not in that format
 * - returns -2 and leaves r unchanged if the denominator is zero
 * - returns 0 otherwise
 */
int rat_set_from_string(rat_t r, const char* s) {
  int retval;
  mpq_t q0;

  mpq_init2(q0, 64);
  // GMP rejects an initial '+' so skip it
  if (*s == '+') s++;
  if (mpq_set_str(q0, s, 10) < 0) {
    retval = -1;
    goto clean_up;
  }
  if (mpz_sgn(mpq_denref(q0)) == 0) {
    retval = -2;  // the denominator is zero
    goto clean_up;
  }
  mpq_canonicalize(q0);
  retval = rat_set_q0(r, &q0);

clean_up:
  mpq_clear(q0);
  return retval;
}

/*
 * Conversion from a string using the given base.
 * Base is interpreted as in gmp: either 0 or an integer from 2 to 36.
 * If base = 0, then the base is determined independently for the
 * numerator and denominator from the first characters.
 * Prefixes  0x or 0b or 0 indicate base 16, 2, or 8,
 * otherwise, the base is 10.
 */
int rat_set_from_string_base(rat_t r, const char* s, int32_t base) {
  int retval;
  mpq_t q0;

  mpq_init2(q0, 64);

  // GMP rejects an initial '+' so skip it
  if (*s == '+') s++;
  assert(0 == base || (2 <= base && base <= 36));
  if (mpq_set_str(q0, s, base) < 0) {
    retval = -1;
    goto clean_up;
  }
  if (mpz_sgn(mpq_denref(q0)) == 0) {
    retval = -2;  // the denominator is zero
    goto clean_up;
  }
  mpq_canonicalize(q0);
  retval = rat_set_q0(r, &q0);

clean_up:

  mpq_clear(q0);
  return retval;
}

/*
 * Portable replacement for strtol that we use below to parse an exponent.
 * The input string must have a non-empty prefix of the form
 *     <optional sign> <digits>
 * without space.
 *
 * Portability issues with strtol
 * On some systems strtol("xxx", ..., 10) sets errno to EINVAL
 * On other systems strtol("xxx", ..., 10) returns 0 and doesn't set errno.
 */
static bool parse_exponent(const char* s, long* result) {
  unsigned long x, y;
  int digit;
  ;
  char c;
  bool ok, positive;

  positive = true;
  ok = false;

  c = *s++;
  if (c == '-') {
    positive = false;
    c = *s++;
  } else if (c == '+') {
    c = *s++;
  }

  x = 0;
  while ('0' <= c && c <= '9') {
    ok = true;
    digit = c - '0';
    y = 10 * x + digit;
    if (y < x) {
      // overflow
      ok = false;
      break;
    }
    x = y;
    c = *s++;
  }

  if (ok) {
    if (positive && x <= (unsigned long)LONG_MAX) {
      *result = (long)x;
      return true;
    }
    if (!positive && x <= (unsigned long)LONG_MIN) {
      *result = (long)(-x);
      return true;
    }
  }

  return false;
}

/*
 * Conversion from a string in a floating point format
 * The expected format is one of
 *   <optional sign> <integer part> . <fractional part>
 *   <optional sign> <integer part> <exp> <optional sign> <integer>
 *   <optional sign> <integer part> . <fractional part> <exp> <optional sign> <integer>
 *
 * Where <optional sign> is + or - and <exp> is either 'e' or 'E'
 *
 * - returns -1 and leaves r unchanged if the string is not in that format
 * - returns 0 otherwise
 */
static int _o_q_set_from_float_string(rat_t r, const char* s) {
  size_t len;
  int frac_len, sign;
  long int exponent;
  char *b, c;
  int retval;
  mpz_t z0;
  mpq_t q0;

  mpz_init2(z0, 64);
  mpq_init2(q0, 64);

  len = strlen(s);
  if (len >= (size_t)UINT32_MAX) {
    // just to be safe if s is really long
    EXIT("string too long in _o_q_set_from_float_string");
  }
  resize_string_buffer(len + 1);
  c = *s++;

  // get sign
  sign = 1;
  if (c == '-') {
    sign = -1;
    c = *s++;
  } else if (c == '+') {
    c = *s++;
  }

  // copy integer part into buffer.
  b = string_buffer;
  while ('0' <= c && c <= '9') {
    *b++ = c;
    c = *s++;
  }

  // copy fractional part and count its length
  frac_len = 0;
  if (c == '.') {
    c = *s++;
    while ('0' <= c && c <= '9') {
      frac_len++;
      *b++ = c;
      c = *s++;
    }
  }
  *b = '\0';  // end of buffer

  // check and read exponent
  exponent = 0;
  if (c == 'e' || c == 'E') {
    if (!parse_exponent(s, &exponent)) {
      retval = -1;
      goto clean_up;
    }
  }

#if 0
  printf("--> Float conversion\n");
  printf("--> sign = %d\n", sign);
  printf("--> mantissa = %s\n", string_buffer);
  printf("--> frac_len = %d\n", frac_len);
  printf("--> exponent = %ld\n", exponent);
#endif

  mpq_set_ui(q0, 0, 1);

  // set numerator
  if (mpz_set_str(mpq_numref(q0), string_buffer, 10) < 0) {
    retval = -1;
    goto clean_up;
  }
  if (sign < 0) {
    mpq_neg(q0, q0);
  }

  // multiply by 10^exponent
  exponent -= frac_len;
  if (exponent > 0) {
    mpz_ui_pow_ui(z0, 10, exponent);
    mpz_mul(mpq_numref(q0), mpq_numref(q0), z0);
  } else if (exponent < 0) {
    // this works even if exponent == LONG_MIN.
    mpz_ui_pow_ui(mpq_denref(q0), 10UL, (unsigned long)(-exponent));
    mpq_canonicalize(q0);
  }

  retval = rat_set_q0(r, &q0);

clean_up:

  mpz_clear(z0);
  mpq_clear(q0);

  return retval;
}

int rat_set_from_float_string(rat_t r, const char* s) { return _o_q_set_from_float_string(r, s); }

/****************
 *  ARITHMETIC  *
 ***************/

void rat_add(rat_t r1, const rat_t r2, const rat_t r3) {
  if (r1 == r2) {
    rat_add_self(r1, r3);
  } else if (r1 == r3) {
    rat_add_self(r1, r2);
  } else {
    rat_set(r1, r2);
    rat_add_self(r1, r3);
  }
}
void rat_sub(rat_t r1, const rat_t r2, const rat_t r3) {
  if (r1 == r2) {
    rat_sub_self(r1, r3);
  } else if (r1 == r3) {
    rat_neg_self(r1);
    rat_add_self(r1, r2);
  } else {
    rat_set(r1, r2);
    rat_sub_self(r1, r3);
  }
}
void rat_mul(rat_t r1, const rat_t r2, const rat_t r3) {
  if (r1 == r2) {
    rat_mul_self(r1, r3);
  } else if (r1 == r3) {
    rat_mul_self(r1, r2);
  } else {
    rat_set(r1, r2);
    rat_mul_self(r1, r3);
  }
}
void rat_div(rat_t r1, const rat_t r2, const rat_t r3) {
  if (r1 == r2) {
    rat_div_self(r1, r3);
  } else if (r1 == r3) {
    rat_inv_self(r1);
    rat_mul_self(r1, r2);
  } else {
    rat_set(r1, r2);
    rat_div_self(r1, r3);
  }
}
void rat_neg(rat_t r1, const rat_t r2) {
  if (r1 == r2) {
    rat_neg_self(r1);
  } else {
    rat_set_neg(r1, r2);
  }
}
void rat_inv(rat_t r1, const rat_t r2) {
  if (r1 == r2) {
    rat_inv_self(r1);
  } else {
    rat_set(r1, r2);
    rat_inv_self(r1);
  }
}

/*
 * Add r2 to r1
 */
void rat_add_self(rat_t r1, const rat_t r2) {
  uint64_t den;
  int64_t num;
  mpq_ptr q1, q2;

  if (r1->s.den == RAT_ONE_DEN && r2->s.den == RAT_ONE_DEN) {
    assert(rat_is_rat32(r2) && rat_is_rat32(r1));
    r1->s.num += r2->s.num;
    if (r1->s.num < MIN_NUMERATOR || r1->s.num > MAX_NUMERATOR) {
      convert_to_gmp(r1);
    }
    return;
  }

  if (rat_is_gmp(r2)) {
    if (!rat_is_gmp(r1)) convert_to_gmp(r1);
    q1 = rat_get_mpq(r1);
    q2 = rat_get_mpq(r2);
    mpq_add(q1, q1, q2);
  } else if (rat_is_gmp(r1)) {
    q1 = rat_get_mpq(r1);
    mpq_add_si(q1, rat_get_num(r2), rat_get_den(r2));
  } else {
    den = rat_get_den(r1) * ((uint64_t)rat_get_den(r2));
    num = rat_get_den(r1) * ((int64_t)rat_get_num(r2)) + rat_get_den(r2) * ((int64_t)rat_get_num(r1));
    rat_set_int64(r1, num, den);
  }
}

/*
 * Subtract r2 from r1
 */
void rat_sub_self(rat_t r1, const rat_t r2) {
  uint64_t den;
  int64_t num;
  mpq_ptr q1, q2;

  if (r1->s.den == RAT_ONE_DEN && r2->s.den == RAT_ONE_DEN) {
    assert(rat_is_rat32(r1) && rat_is_rat32(r2));
    r1->s.num -= r2->s.num;
    if (r1->s.num < MIN_NUMERATOR || r1->s.num > MAX_NUMERATOR) {
      convert_to_gmp(r1);
    }
    return;
  }

  if (rat_is_gmp(r2)) {
    if (!rat_is_gmp(r1)) convert_to_gmp(r1);
    q1 = rat_get_mpq(r1);
    q2 = rat_get_mpq(r2);
    mpq_sub(q1, q1, q2);
  } else if (rat_is_gmp(r1)) {
    q1 = rat_get_mpq(r1);
    mpq_sub_si(q1, rat_get_num(r2), rat_get_den(r2));
  } else {
    den = rat_get_den(r1) * ((uint64_t)rat_get_den(r2));
    num = rat_get_den(r2) * ((int64_t)rat_get_num(r1)) - rat_get_den(r1) * ((int64_t)rat_get_num(r2));
    rat_set_int64(r1, num, den);
  }
}

/*
 * Negate r
 */
void rat_neg_self(rat_t r) {
  if (rat_is_gmp(r)) {
    mpq_ptr q = rat_get_mpq(r);
    mpq_neg(q, q);
  } else {
    r->s.num = -r->s.num;
  }
}

/*
 * Invert r
 */
void rat_inv_self(rat_t r) {
  uint32_t abs_num;

  if (rat_is_gmp(r)) {
    mpq_ptr q = rat_get_mpq(r);
    mpq_inv(q, q);
  } else if (r->s.num < 0) {
    abs_num = (uint32_t)-r->s.num;
    set_rat32(r, -rat_get_den(r), abs_num);
  } else if (r->s.num > 0) {
    abs_num = (uint32_t)r->s.num;
    set_rat32(r, rat_get_den(r), abs_num);
  } else {
    division_by_zero();
  }
}

/*
 * Invert r mod m
 */
void rat_inv_mod(rat_t r, const rat_t mod) {
  assert(rat_is_integer(r) && rat_is_integer(mod) && rat_is_pos(mod));

  rat_s tmp, *mm;
  if (rat_is_rat32(mod)) {
    tmp = *mod;
    mm = &tmp;
    rat_denormalize(mm);
  } else {
    mm = (rat_ptr)mod;
  }

  rat_denormalize(r);
  assert(rat_is_gmp(r) && rat_is_gmp(mm));
  mpz_ptr z = rat_get_mpq_num(r);
  mpz_ptr m = rat_get_mpq_num(mm);
  mpz_invert(z, z, m);
  rat_normalize(r);

  if (mm != mod) rat_clear(mm);
}

void rat_add_si(rat_t r1, const int64_t num, const uint64_t den) {
  if (rat_is_gmp(r1)) {
    mpq_ptr q1 = rat_get_mpq(r1);
    mpq_add_si(q1, num, den);
  } else {
    rat_t tmp;
    rat_init(tmp);
    rat_set_int64(tmp, num, den);
    rat_add_self(r1, tmp);
    rat_clear(tmp);
  }
}
void rat_sub_si(rat_t r1, const int64_t num, const uint64_t den) {
  if (rat_is_gmp(r1)) {
    mpq_ptr q1 = rat_get_mpq(r1);
    mpq_sub_si(q1, num, den);
  } else {
    rat_t tmp;
    rat_init(tmp);
    rat_set_int64(tmp, num, den);
    rat_sub_self(r1, tmp);
    rat_clear(tmp);
  }
}
void rat_mul_si(rat_t r1, const int64_t num, const uint64_t den) {
  if (rat_is_gmp(r1)) {
    mpq_ptr q1 = rat_get_mpq(r1);
    mpq_mul_si(q1, num, den);
  } else {
    rat_t tmp;
    rat_init(tmp);
    rat_set_int64(tmp, num, den);
    rat_mul_self(r1, tmp);
    rat_clear(tmp);
  }
}
void rat_div_si(rat_t r1, const int64_t num, const uint64_t den) {
  if (rat_is_gmp(r1)) {
    mpq_ptr q1 = rat_get_mpq(r1);
    mpq_div_si(q1, num, den);
  } else {
    rat_t tmp;
    rat_init(tmp);
    rat_set_int64(tmp, num, den);
    rat_div_self(r1, tmp);
    rat_clear(tmp);
  }
}

extern void rat_add_mpq(rat_t r1, const mpq_t q) {
  if (rat_is_gmp(r1)) {
    mpq_ptr q1 = rat_get_mpq(r1);
    mpq_add(q1, q1, q);
  } else {
    rat_t tmp;
    rat_init(tmp);
    rat_set_mpq(tmp, q);
    rat_add_self(r1, tmp);
    rat_clear(tmp);
  }
}
extern void rat_sub_mpq(rat_t r1, const mpq_t q) {
  if (rat_is_gmp(r1)) {
    mpq_ptr q1 = rat_get_mpq(r1);
    mpq_sub(q1, q1, q);
  } else {
    rat_t tmp;
    rat_init(tmp);
    rat_set_mpq(tmp, q);
    rat_sub_self(r1, tmp);
    rat_clear(tmp);
  }
}
extern void rat_mul_mpq(rat_t r1, const mpq_t q) {
  if (rat_is_gmp(r1)) {
    mpq_ptr q1 = rat_get_mpq(r1);
    mpq_mul(q1, q1, q);
  } else {
    rat_t tmp;
    rat_init(tmp);
    rat_set_mpq(tmp, q);
    rat_mul_self(r1, tmp);
    rat_clear(tmp);
  }
}
extern void rat_div_mpq(rat_t r1, const mpq_t q) {
  if (rat_is_gmp(r1)) {
    mpq_ptr q1 = rat_get_mpq(r1);
    mpq_div(q1, q1, q);
  } else {
    rat_t tmp;
    rat_init(tmp);
    rat_set_mpq(tmp, q);
    rat_div_self(r1, tmp);
    rat_clear(tmp);
  }
}

/*
 * Multiply r1 by r2
 */
void rat_mul_self(rat_t r1, const rat_t r2) {
  uint64_t den;
  int64_t num;
  mpq_ptr q1, q2;

  if (r1->s.den == RAT_ONE_DEN && r2->s.den == RAT_ONE_DEN) {
    assert(rat_is_rat32(r1) && rat_is_rat32(r2));
    num = r1->s.num * ((int64_t)r2->s.num);
    if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR) {
      r1->s.num = (int32_t)num;
    } else {
      set_to_gmp64(r1, num);
    }
    return;
  }

  if (rat_is_gmp(r2)) {
    if (rat_is_rat32(r1)) {
      convert_to_gmp(r1);
    }
    q1 = rat_get_mpq(r1);
    q2 = rat_get_mpq(r2);
    mpq_mul(q1, q1, q2);
  } else if (rat_is_gmp(r1)) {
    q1 = rat_get_mpq(r1);
    mpq_mul_si(q1, rat_get_num(r2), rat_get_den(r2));
  } else {
    den = rat_get_den(r1) * ((uint64_t)rat_get_den(r2));
    num = rat_get_num(r1) * ((int64_t)rat_get_num(r2));
    rat_set_int64(r1, num, den);
  }
}

/*
 * Divide r1 by r2
 */
void rat_div_self(rat_t r1, const rat_t r2) {
  uint64_t den;
  int64_t num;
  mpq_ptr q1, q2;

  if (rat_is_gmp(r2)) {
    if (rat_is_rat32(r1)) convert_to_gmp(r1);
    q1 = rat_get_mpq(r1);
    q2 = rat_get_mpq(r2);
    mpq_div(q1, q1, q2);
  } else if (rat_is_gmp(r1)) {
    if (rat_get_num(r2) == 0) {
      division_by_zero();
    } else {
      q1 = rat_get_mpq(r1);
      mpq_div_si(q1, rat_get_num(r2), rat_get_den(r2));
    }

  } else if (rat_get_num(r2) > 0) {
    den = rat_get_den(r1) * ((uint64_t)rat_get_num(r2));
    num = rat_get_num(r1) * ((int64_t)rat_get_den(r2));
    rat_set_int64(r1, num, den);

  } else if (rat_get_num(r2) < 0) {
    den = rat_get_den(r1) * ((uint64_t)(-rat_get_num(r2)));
    num = rat_get_num(r1) * (-((int64_t)rat_get_den(r2)));
    rat_set_int64(r1, num, den);

  } else {
    division_by_zero();
  }
}

/*
 * Add r2 * r3 to  r1
 */
void rat_addmul(rat_t r1, const rat_t r2, const rat_t r3) {
  int64_t num;
  rat_t tmp;

  if (r1->s.den == RAT_ONE_DEN && r2->s.den == RAT_ONE_DEN && r3->s.den == RAT_ONE_DEN) {
    assert(rat_is_rat32(r1) && rat_is_rat32(r2) && rat_is_rat32(r3));

    num = rat_get_num(r1) + rat_get_num(r2) * ((int64_t)rat_get_num(r3));
    if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR) {
      r1->s.num = (int32_t)num;
    } else {
      set_to_gmp64(r1, num);
    }
    return;
  }

  rat_init(tmp);
  rat_set(tmp, r2);
  rat_mul_self(tmp, r3);
  rat_add_self(r1, tmp);
  rat_clear(tmp);
}

/*
 * Subtract r2 * r3 from r1
 */
void rat_submul(rat_t r1, const rat_t r2, const rat_t r3) {
  int64_t num;
  rat_t tmp;

  if (r1->s.den == RAT_ONE_DEN && r2->s.den == RAT_ONE_DEN && r3->s.den == RAT_ONE_DEN) {
    assert(rat_is_rat32(r1) && rat_is_rat32(r2) && rat_is_rat32(r3));

    num = rat_get_num(r1) - rat_get_num(r2) * ((int64_t)rat_get_num(r3));
    if (MIN_NUMERATOR <= num && num <= MAX_NUMERATOR) {
      r1->s.num = (int32_t)num;
    } else {
      set_to_gmp64(r1, num);
    }
    return;
  }

  rat_init(tmp);
  rat_set(tmp, r2);
  rat_mul_self(tmp, r3);
  rat_sub_self(r1, tmp);
  rat_clear(tmp);
}

/*
 * Increment: add one to r1
 */
void rat_add_one(rat_t r1) {
  mpq_ptr q;

  if (rat_is_gmp(r1)) {
    q = rat_get_mpq(r1);
    mpz_add(mpq_numref(q), mpq_numref(q), mpq_denref(q));
  } else {
    r1->s.num += rat_get_den(r1);
    if (r1->s.num > MAX_NUMERATOR) {
      convert_to_gmp(r1);
    }
  }
}

/*
 * Decrement: subtract one from r1
 */
void rat_sub_one(rat_t r1) {
  mpq_ptr q;

  if (rat_is_gmp(r1)) {
    q = rat_get_mpq(r1);
    mpz_sub(mpq_numref(q), mpq_numref(q), mpq_denref(q));
  } else {
    r1->s.num -= rat_get_den(r1);
    if (r1->s.num < MIN_NUMERATOR) {
      convert_to_gmp(r1);
    }
  }
}

/***********************
 *  CEILING AND FLOOR  *
 **********************/

// set r to floor(r);
void rat_floor(rat_t r1, const rat_t r2) {
  if (r1 == r2) {
    rat_floor_self(r1);
  } else {
    rat_set(r1, r2);
    rat_floor_self(r1);
  }
}

void rat_ceil(rat_t r1, const rat_t r2) {
  if (r1 == r2) {
    rat_ceil_self(r1);
  } else {
    rat_set(r1, r2);
    rat_ceil_self(r1);
  }
}

void rat_floor_self(rat_t r) {
  int32_t n;

  if (rat_is_integer(r)) return;

  if (rat_is_gmp(r)) {
    mpq_ptr q = rat_get_mpq(r);
    mpz_fdiv_q(mpq_numref(q), mpq_numref(q), mpq_denref(q));
    mpz_set_ui(mpq_denref(q), 1UL);
  } else {
    n = r->s.num / (int32_t)rat_get_den(r);
    if (r->s.num < 0) n--;
    set_rat32(r, n, 1);
  }
}

// set r to ceil(r)
void rat_ceil_self(rat_t r) {
  int32_t n;

  if (rat_is_integer(r)) return;

  if (rat_is_gmp(r)) {
    mpq_ptr q = rat_get_mpq(r);
    mpz_cdiv_q(mpq_numref(q), mpq_numref(q), mpq_denref(q));
    mpz_set_ui(mpq_denref(q), 1UL);
  } else {
    n = r->s.num / (int32_t)rat_get_den(r);
    if (r->s.num > 0) n++;
    set_rat32(r, n, 1);
  }
}

/*******************
 *  EXPONENTATION  *
 ******************/

/*
 * Store r1 * (r2 ^ n) into r1
 */
void rat_mulexp(rat_t r1, const rat_t r2, uint32_t n) {
  rat_t aux;

  if (n <= 3) {
    // small exponent:
    switch (n) {
      case 3:
        rat_mul_self(r1, r2);
      case 2:
        rat_mul_self(r1, r2);
      case 1:
        rat_mul_self(r1, r2);
      case 0:
        break;  // do nothing
    }

  } else {
    rat_init(aux);
    rat_set(aux, r2);

    // compute r1 * aux^n
    for (;;) {
      assert(n > 0);
      if ((n & 1) != 0) {
        rat_mul_self(r1, aux);
      }
      n >>= 1;
      if (n == 0) break;
      rat_mul_self(aux, aux);  // this should work
    }

    rat_clear(aux);
  }
}

/*****************
 *  LCM and GCD  *
 ****************/

/*
 * Get the absolute value of x (converted to uint32_t)
 */
static inline uint32_t abs32(int32_t x) { return (x >= 0) ? x : -x; }

/*
 * Store lcm(r1, r2) into r1
 * - r1 and r2 must be integer
 * - the result is always positive
 */
void rat_lcm(rat_t r1, const rat_t r2) {
  uint32_t a, b;
  uint64_t d;
  mpq_ptr q1, q2;

  if (rat_is_rat32(r2)) {
    if (rat_is_rat32(r1)) {
      // both r1 and r2 are 32bit integers
      a = abs32(r1->s.num);
      b = abs32(r2->s.num);
      d = ((uint64_t)a) * ((uint64_t)(b / gcd32(a, b)));
      if (d <= MAX_NUMERATOR) {
        set_rat32(r1, d, 1);
      } else {
        set_to_gmp64(r1, d);
      }
    } else {
      // r2 is 32bits, r1 is gmp
      b = abs32(r2->s.num);
      q1 = rat_get_mpq(r1);
      mpz_lcm_ui(mpq_numref(q1), mpq_numref(q1), b);
    }
  } else {
    // r2 is a gmp rational
    if (rat_is_rat32(r1)) convert_to_gmp(r1);
    q1 = rat_get_mpq(r1);
    q2 = rat_get_mpq(r2);
    mpz_lcm(mpq_numref(q1), mpq_numref(q1), mpq_numref(q2));
  }

  assert(rat_is_pos(r1));
}

/*
 * Store gcd(r1, r2) into r1
 * - r1 and r2 must be integer and non-zero
 * - the result is positive
 */
void rat_gcd(rat_t r1, const rat_t r2) {
  uint32_t a, b, d;
  mpq_ptr q1, q2;

  if (rat_is_rat32(r2)) {
    if (rat_is_rat32(r1)) {
      // r1 and r2 are small integers
      a = abs32(r1->s.num);
      b = abs32(r2->s.num);
      d = gcd32(a, b);
    } else {
      // r1 is gmp, r2 is a small integer
      b = abs32(r2->s.num);
      q1 = rat_get_mpq(r1);
      d = mpz_gcd_ui(NULL, mpq_numref(q1), b);
      release_mpq(r1);
    }
    assert(d <= MAX_NUMERATOR);
    set_rat32(r1, d, 1);
  } else {
    if (rat_is_rat32(r1)) {
      // r1 is a small integer, r2 is a gmp number
      a = abs32(r1->s.num);
      q2 = rat_get_mpq(r2);
      d = mpz_gcd_ui(NULL, mpq_numref(q2), a);
      assert(d <= MAX_NUMERATOR);
      set_rat32(r1, d, 1);
    } else {
      // both are gmp numbers
      q1 = rat_get_mpq(r1);
      q2 = rat_get_mpq(r2);
      mpz_gcd(mpq_numref(q1), mpq_numref(q1), mpq_numref(q2));
    }
  }
}

/************************
 *  EUCLIDEAN DIVISION  *
 ***********************/

/*
 * Quotient/remainder in Euclidean division
 * - r1 and r2 must be integer
 * - r2 must be positive
 */
void rat_integer_div(rat_t r1, const rat_t r2) {
  int32_t n;
  mpq_ptr q1, q2;

  if (rat_is_rat32(r2)) {
    assert(r2->s.den == RAT_ONE_DEN && r2->s.num > 0);
    if (rat_is_rat32(r1)) {
      // both r1 and r2 are small integers
      n = r1->s.num % r2->s.num;  // remainder: n has the same sign as r1 (or n == 0)
      r1->s.num /= r2->s.num;     // quotient: same sign as r1, rounded towards 0
      if (n < 0) {
        r1->s.num--;
      }
    } else {
      // r1 is gmp, r2 is a small integer
      q1 = rat_get_mpq(r1);
      mpz_fdiv_q_ui(mpq_numref(q1), mpq_numref(q1), r2->s.num);
      assert(mpq_is_integer(q1));
    }
  } else {
    q2 = rat_get_mpq(r2);
    assert(mpq_is_integer(q2) && mpq_sgn(q2) > 0);
    if (rat_is_rat32(r1)) {
      /*
       * r1 is a small integer, r2 is a gmp rational
       * since r2 is normalized and positive, we have r2 > abs(r1)
       * so r1 div r2 = 0 if r1 >= 0 or -1 if r1 < 0
       */
      assert(r1->s.den == RAT_ONE_DEN);
      if (r1->s.num >= 0) {
        r1->s.num = 0;
      } else {
        r1->s.num = -1;
      }
    } else {
      // both r1 and r2 are gmp rationals
      q1 = rat_get_mpq(r1);
      mpz_fdiv_q(mpq_numref(q1), mpq_numref(q1), mpq_numref(q2));
      assert(mpq_is_integer(q1));
    }
  }
}

/*
 * Assign the remainder of r1 divided by r2 to r1
 * - both r1 and r2 must be integer
 * - r2 must be positive
 */
void rat_integer_rem(rat_t r1, const rat_t r2) {
  int32_t n;
  mpq_ptr q1, q2;

  if (rat_is_rat32(r2)) {
    assert(r2->s.den == RAT_ONE_DEN && r2->s.num > 0);
    if (rat_is_rat32(r1)) {
      /*
       * both r1 and r2 are small integers
       * Note: the result of (r1->num % r2->num) has the same sign as r1->num
       */
      n = r1->s.num % r2->s.num;
      if (n < 0) {
        n += r2->s.num;
      }
      assert(0 <= n && n < r2->s.num);
      r1->s.num = n;
    } else {
      // r1 is gmp, r2 is a small integer
      q1 = rat_get_mpq(r1);
      n = mpz_fdiv_ui(mpq_numref(q1), r2->s.num);
      assert(0 <= n && n <= MAX_NUMERATOR);
      release_mpq(r1);
      set_rat32(r1, n, 1);
    }
  } else {
    q2 = rat_get_mpq(r2);

    assert(mpq_is_integer(q2) && mpq_sgn(q2) > 0);
    if (rat_is_rat32(r1)) {
      /*
       * r1 is a small integer, r2 is a gmp rational
       * since r2 is normalized and positive, we have r2 > abs(r1)
       * so r1 mod r2 = r1 if r1 >= 0
       * or r1 mod r2 = (r2 + r1) if r1 < 0
       */
      assert(r1->s.den == RAT_ONE_DEN);
      if (r1->s.num < 0) {
        mpq_ptr q = new_mpq();
        mpq_set_si(q, r1->s.num, 1UL);
        mpz_add(mpq_numref(q), mpq_numref(q), mpq_numref(q2));
        rat_set_gmp(r1, q);
        assert(mpq_is_integer(q) && mpq_sgn(q) > 0);
      }

    } else {
      // both r1 and r2 are gmp rationals
      q1 = rat_get_mpq(r1);
      mpz_fdiv_r(mpq_numref(q1), mpq_numref(q1), mpq_numref(q2));
      assert(mpq_is_integer(q1));
    }
  }
}

/*
 * Check whether r1 divides r2
 *
 * Both r1 and r2 must be integers and r1 must be non-zero
 */
bool rat_integer_divides(rat_t r1, const rat_t r2) {
  uint32_t aux;
  mpq_ptr q1, q2;

  rat_normalize(r1);

  if (rat_is_gmp(r1)) {
    if (rat_is_gmp(r2)) {
      q1 = rat_get_mpq(r1);
      q2 = rat_get_mpq(r2);
      return mpz_divisible_p(mpq_numref(q2), mpq_numref(q1));
    } else {
      return false;  // abs(r1) > abs(r2) so r1 can't divide r2
    }

  } else {
    assert(r1->s.den == RAT_ONE_DEN);
    aux = abs32(r1->s.num);
    if (rat_is_gmp(r2)) {
      q2 = rat_get_mpq(r2);
      return mpz_divisible_ui_p(mpq_numref(q2), aux);
    } else {
      return abs32(r2->s.num) % aux == 0;
    }
  }
}

/*
 * Check whether r2/r1 is an integer
 * - r1 must be non-zero
 */
bool rat_divides(const rat_t r1, const rat_t r2) {
  rat_t aux;
  bool divides;

  if (r1->s.den == RAT_ONE_DEN && (r1->s.num == 1 || r1->s.num == -1)) {
    assert(rat_is_rat32(r1));
    // r1 is +1 or -1
    return true;
  }

  rat_init(aux);
  rat_set(aux, r2);
  rat_div_self(aux, r1);
  divides = rat_is_integer(aux);
  rat_clear(aux);

  return divides;
}

/***********************************
 *  SMT2 version of (divides x y)  *
 **********************************/

/*
 * (divides r1 r2) is (exist (n::int) r2 = n * r1)
 *
 * This is the same as rat_divides(r1, r2) except that the
 * definition allows r1 to be zero. In this case,
 *  (divides 0 r2) is (r2 == 0)
 */
bool rat_smt2_divides(const rat_t r1, const rat_t r2) {
  if (rat_is_zero(r1)) {
    return rat_is_zero(r2);
  } else {
    return rat_divides(r1, r2);
  }
}

/**********************
 *  GENERALIZED LCM   *
 *********************/

/*
 * This computes the LCM of r1 and r2 for arbitrary rationals:
 * - if r1 is (a1/b1) and r2 is (a2/b2) then the result is
 *    lcm(a1, a2)/gcd(b1, b2).
 * - the result is stored in r1
 */
void rat_generalized_lcm(rat_t r1, const rat_t r2) {
  rat_t a1, b1;
  rat_t a2, b2;

  if (rat_is_integer(r1) && rat_is_integer(r2)) {
    rat_lcm(r1, r2);
  } else {
    rat_init(a1);
    rat_set_num(a1, r1);
    rat_init(b1);
    rat_set_den(b1, r1);

    rat_init(a2);
    rat_set_num(a2, r2);
    rat_init(b2);
    rat_set_den(b2, r2);

    rat_lcm(a1, a2);  // a1 := lcm(a1, a2)
    rat_gcd(b1, b2);  // b1 := gcd(b1, b2)

    // copy the result in r1
    rat_set(r1, a1);
    rat_div_self(r1, b1);

    rat_clear(a1);
    rat_clear(b1);
    rat_clear(a2);
    rat_clear(b2);
  }
}

/*
 * This computes the GCD of r1 and r2 for arbitrary non-zero rationals:
 * - if r1 is (a1/b1) and r2 is (a2/b2) then the result is
 *    gcd(a1, a2)/lcm(b1, b2).
 * - the result is stored in r1
 */
void rat_generalized_gcd(rat_t r1, const rat_t r2) {
  rat_t a1, b1;
  rat_t a2, b2;

  if (rat_is_integer(r1) && rat_is_integer(r2)) {
    rat_lcm(r1, r2);
  } else {
    rat_init(a1);
    rat_set_num(a1, r1);
    rat_init(b1);
    rat_set_den(b1, r1);

    rat_init(a2);
    rat_set_num(a2, r2);
    rat_init(b2);
    rat_set_den(b2, r2);

    rat_gcd(a1, a2);  // a1 := gcd(a1, a2)
    rat_lcm(b1, b2);  // b1 := lcm(b1, b2)

    // copy the result in r1
    rat_set(r1, a1);
    rat_div_self(r1, b1);

    rat_clear(a1);
    rat_clear(b1);
    rat_clear(a2);
    rat_clear(b2);
  }
}

/**********************
 *  SMT2 DIV AND MOD  *
 *********************/

/*
 * SMT-LIB 2.0 definitions for div and mod:
 * - if y > 0 then div(x, y) is floor(x/y)
 * - if y < 0 then div(x, y) is ceil(x/y)
 * - 0 <= mod(x, y) < y
 * - x = y * div(x, y) + mod(x, y)
 * These operations are defined for any x and non-zero y.
 * The terms x and y are not required to be integers.
 *
 * - rat_smt2_div(q, x, y) stores (div x y) in q
 * - rat_smt2_mod(q, x, y) stores (mod x y) in q
 *
 * For both functions, y must not be zero.
 */
void rat_smt2_div(rat_t q, const rat_t x, const rat_t y) {
  assert(rat_is_nonzero(y));

  rat_set(q, x);
  rat_div_self(q, y);  // q := x/y
  if (rat_is_pos(y)) {
    rat_floor_self(q);  // round down
  } else {
    rat_ceil_self(q);  // round up
  }
}

/*
 * For debugging: check that 0 <= r < abs(y)
 */
#ifndef NDEBUG
static bool plausible_mod(const rat_t r, const rat_t y) {
  rat_t aux;
  bool ok;

  assert(rat_is_nonzero(y));

  rat_init(aux);
  if (rat_is_pos(y)) {
    rat_set(aux, y);
  } else {
    rat_set_neg(aux, y);
  }
  rat_normalize(aux);

  ok = rat_is_nonneg(r) && rat_lt(r, aux);

  rat_clear(aux);

  return ok;
}
#endif

void rat_smt2_mod(rat_t q, const rat_t x, const rat_t y) {
  assert(rat_is_nonzero(y));

  rat_smt2_div(q, x, y);  // q := (div x y)
  rat_mul_self(q, y);     // q := y * (div x y)
  rat_sub_self(q, x);     // q := - x + y * (div x y)
  rat_neg_self(q);        // q := x - y * (div x y) = (mod x y)

  assert(plausible_mod(q, y));
}

/*****************
 *  COMPARISONS  *
 ****************/

/*
 * Compare r1 and r2
 * - returns a negative number if r1 < r2
 * - returns 0 if r1 = r2
 * - returns a positive number if r1 > r2
 */
int rat_cmp(const rat_t r1, const rat_t r2) {
  int64_t num;
  mpq_ptr q1, q2;

  if (r1->s.den == RAT_ONE_DEN && r2->s.den == RAT_ONE_DEN) {
    assert(rat_is_rat32(r1) && rat_is_rat32(r2));
    return r1->s.num - r2->s.num;
  }

  if (rat_is_gmp(r1)) {
    if (rat_is_gmp(r2)) {
      q1 = rat_get_mpq(r1);
      q2 = rat_get_mpq(r2);
      return mpq_cmp(q1, q2);
    } else {
      q1 = rat_get_mpq(r1);
      return mpq_cmp_si(q1, rat_get_num(r2), rat_get_den(r2));
    }
  } else {
    if (rat_is_gmp(r2)) {
      q2 = rat_get_mpq(r2);
      return -mpq_cmp_si(q2, rat_get_num(r1), rat_get_den(r1));
    } else {
      num = rat_get_den(r2) * ((int64_t)rat_get_num(r1)) - rat_get_den(r1) * ((int64_t)rat_get_num(r2));
      return (num < 0 ? -1 : (num > 0));
    }
  }
}

/*
 * Compare r1 and num/den
 */
int rat_cmp_int32(const rat_t r1, int32_t num, uint32_t den) {
  int64_t nn;
  mpq_ptr q1;

  if (rat_is_gmp(r1)) {
    q1 = rat_get_mpq(r1);
    return mpq_cmp_si(q1, num, den);
  } else {
    nn = den * ((int64_t)r1->s.num) - rat_get_den(r1) * ((int64_t)num);
    return (nn < 0 ? -1 : (nn > 0));
  }
}

int rat_cmp_int64(const rat_t r1, int64_t num, uint64_t den) {
  int retval;
  mpq_ptr q1;
  mpq_t q0;

  mpq_init2(q0, 64);
  mpq_set_int64(q0, num, den);
  mpq_canonicalize(q0);
  if (rat_is_gmp(r1)) {
    q1 = rat_get_mpq(r1);
    retval = mpq_cmp(q1, q0);
  } else {
    retval = -mpq_cmp_si(q0, r1->s.num, rat_get_den(r1));
  }
  mpq_clear(q0);
  return retval;
}

/*
 * Check whether r1 and r2 are opposite
 */
bool rat_opposite(const rat_t r1, const rat_t r2) {
  rat_t aux;
  bool result;

  if (r1->s.den == RAT_ONE_DEN && r2->s.den == RAT_ONE_DEN) {
    assert(rat_is_rat32(r1) && rat_is_rat32(r2));
    return r1->s.num + r2->s.num == 0;
  }

  rat_init(aux);
  rat_set(aux, r1);
  rat_add_self(aux, r2);
  result = rat_is_zero(aux);
  rat_clear(aux);

  return result;
}

/***********************************************
 *  CONVERSIONS FROM RATIONALS TO OTHER TYPES  *
 **********************************************/

/*
 * Convert r to a 32bit signed integer
 * - return false if r is not an integer or does not fit in 32 bits
 */
bool rat_get32(rat_t r, int32_t* v) {
  uint32_t d;
  mpq_ptr q;

  if (r->s.den == RAT_ONE_DEN) {
    assert(rat_is_rat32(r));
    *v = r->s.num;
    return true;
  } else if (rat_is_gmp(r)) {
    q = rat_get_mpq(r);
    if (mpq_fits_int32(q)) {
      mpq_get_int32(q, v, &d);
      return d == 1;
    }
  }
  return false;
}

/*
 * Convert r to a 64bit signed integer v
 * - return false if r is not an integer or does not fit in 64bits
 */
bool rat_get64(rat_t r, int64_t* v) {
  uint64_t d;
  mpq_ptr q;

  if (r->s.den == RAT_ONE_DEN) {
    assert(rat_is_rat32(r));
    *v = r->s.num;
    return true;
  } else if (rat_is_gmp(r)) {
    q = rat_get_mpq(r);
    if (mpq_fits_int64(q)) {
      mpq_get_int64(q, v, &d);
      return d == 1;
    }
  }
  return false;
}

/*
 * Convert r to a pair of 32 bit integers num/den
 * - return false if the numerator or denominator doesn't fit in 32bits
 */
bool rat_get_int32(rat_t r, int32_t* num, uint32_t* den) {
  mpq_ptr q;

  if (rat_is_rat32(r)) {
    *num = rat_get_num(r);
    *den = rat_get_den(r);
    return true;
  } else {
    q = rat_get_mpq(r);
    if (mpq_fits_int32(q)) {
      mpq_get_int32(q, num, den);
      return true;
    }
  }
  return false;
}

/*
 * Convert r to a pair of 64bit integers num/den
 * - return false if the numerator or denominator doesn't fit in 64bits
 */
bool rat_get_int64(rat_t r, int64_t* num, uint64_t* den) {
  mpq_ptr q;

  if (rat_is_rat32(r)) {
    *num = rat_get_num(r);
    *den = rat_get_den(r);
    return true;
  } else {
    q = rat_get_mpq(r);
    if (mpq_fits_int64(q)) {
      mpq_get_int64(q, num, den);
      return true;
    }
  }
  return false;
}

/*
 * Check whether r can be converted to a 32bit integer,
 * a 64bit integer, or two a pair num/den of 32bit or 64bit integers.
 */

bool rat_is_int32(const rat_t r) {
  return (rat_is_rat32(r) && r->s.den == RAT_ONE_DEN) || (rat_is_gmp(r) && mpq_is_int32(rat_get_mpq(r)));
}

bool rat_is_int64(const rat_t r) {
  return (rat_is_rat32(r) && r->s.den == RAT_ONE_DEN) || (rat_is_gmp(r) && mpq_is_int64(rat_get_mpq(r)));
}

bool rat_fits_int32(const rat_t r) { return rat_is_rat32(r) || mpq_fits_int32(rat_get_mpq(r)); }
bool rat_fits_int64(const rat_t r) { return rat_is_rat32(r) || mpq_fits_int64(rat_get_mpq(r)); }

/*
 * Size estimate
 * - r must be an integer
 * - this returns approximately the number of bits to represent r
 */
uint32_t rat_size(const rat_t r) {
  size_t n;

  n = 32;
  if (rat_is_gmp(r)) {
    n = mpz_size(mpq_numref(rat_get_mpq(r))) * mp_bits_per_limb;
    if (n > (size_t)UINT32_MAX) {
      n = UINT32_MAX;
    }
  }
  return (uint32_t)n;
}

/*
 * Convert r to a GMP integer
 * - return false if r is not an integer
 */
bool rat_to_mpz(const rat_t r, mpz_t z) {
  if (r->s.den == RAT_ONE_DEN) {
    assert(rat_is_rat32(r));
    mpz_set_si(z, r->s.num);
    return true;
  } else if (rat_is_gmp(r)) {
    mpq_ptr q = rat_get_mpq(r);
    if (mpq_is_integer(q)) {
      mpz_set(z, mpq_numref(q));
      return true;
    }
  }
  return false;
}

/*
 * Convert r to a GMP rational
 */
void rat_to_mpq(const rat_t r, mpq_t q) {
  if (rat_is_gmp(r)) {
    mpq_set(q, rat_get_mpq(r));
  } else {
    mpq_set_int32(q, r->s.num, rat_get_den(r));
  }
}

/*
 * Convert to a double
 */
double rat_to_d(const rat_t r) {
  double retval;
  mpq_t q0;

  mpq_init2(q0, 64);
  rat_to_mpq(r, q0);
  retval = mpq_get_d(q0);
  mpq_clear(q0);

  return retval;
}

/**************
 *  PRINTING  *
 *************/

/*
 * Print r
 */
void rat_print(FILE* f, const rat_t r) {
  if (rat_is_gmp(r)) {
    mpq_out_str(f, 10, rat_get_mpq(r));
  } else if (r->s.den != RAT_ONE_DEN) {
    fprintf(f, "%" PRId32 "/%" PRIu32, r->s.num, rat_get_den(r));
  } else {
    fprintf(f, "%" PRId32, r->s.num);
  }
}

/*
 * Print r's absolute value
 */
void rat_print_abs(FILE* f, const rat_t r) {
  mpq_ptr q;
  int32_t abs_num;

  if (rat_is_gmp(r)) {
    q = rat_get_mpq(r);
    if (mpq_sgn(q) < 0) {
      mpq_neg(q, q);
      mpq_out_str(f, 10, q);
      mpq_neg(q, q);
    } else {
      mpq_out_str(f, 10, q);
    }
  } else {
    abs_num = r->s.num;
    if (abs_num < 0) abs_num = -abs_num;

    if (r->s.den != RAT_ONE_DEN) {
      fprintf(f, "%" PRId32 "/%" PRIu32, abs_num, rat_get_den(r));
    } else {
      fprintf(f, "%" PRId32, abs_num);
    }
  }
}

char* rat_get_str(const rat_t r) {
  char* s;

  if (rat_is_gmp(r)) {
    s = mpq_get_str(NULL, 10, rat_get_mpq(r));
  } else if (r->s.den != RAT_ONE_DEN) {
    s = (char*)EGmalloc(32 + 32 + 2);
    sprintf(s, "%" PRId32 "/%" PRIu32, rat_get_num(r), rat_get_den(r));
  } else {
    s = (char*)EGmalloc(32 + 1);
    sprintf(s, "%" PRId32, rat_get_num(r));
  }

  return s;
}

/****************
 *  HASH CODES  *
 ***************/

/* largest prime less than 2^32 */
#define HASH_MODULUS 4294967291UL

/*
 *  we have  HASH_MODULUS > - MIN_NUMERATOR
 *   and HASH_MODULUS > MAX_NUMERATOR
 *
 * if MIN_NUMERATOR <= r->num <= MAX_NUMERATOR then
 * 1) if r->num >= 0 then r_num mod HASH_MODULUS = r_nun
 * 2) if r->num < 0 then r_num mod HASH_MODULUS = r_num + HASH_MODULUS
 *  so  ((uint32_t) r->num) + HASH_MODULUS = (2^32 + r->num) + HASH_MODULUS
 * gives the correct result.
 */
uint32_t rat_hash_numerator(const rat_t r) {
  if (rat_is_gmp(r)) {
    return (uint32_t)mpz_fdiv_ui(mpq_numref(rat_get_mpq(r)), HASH_MODULUS);
  } else if (r->s.num >= 0) {
    return (uint32_t)r->s.num;
  } else {
    return ((uint32_t)r->s.num) + ((uint32_t)HASH_MODULUS);
  }
}

uint32_t rat_hash_denominator(const rat_t r) {
  if (rat_is_gmp(r)) {
    return (uint32_t)mpz_fdiv_ui(mpq_denref(rat_get_mpq(r)), HASH_MODULUS);
  }
  return rat_get_den(r);
}

void rat_hash_decompose(const rat_t r, uint32_t* h_num, uint32_t* h_den) {
  if (rat_is_gmp(r)) {
    mpq_ptr q = rat_get_mpq(r);
    *h_num = (uint32_t)mpz_fdiv_ui(mpq_numref(q), HASH_MODULUS);
    *h_den = (uint32_t)mpz_fdiv_ui(mpq_denref(q), HASH_MODULUS);
  } else if (r->s.num >= 0) {
    *h_num = (uint32_t)r->s.num;
    *h_den = rat_get_den(r);
  } else {
    *h_num = ((uint32_t)r->s.num) + ((uint32_t)HASH_MODULUS);
    *h_den = rat_get_den(r);
  }
}

/***********************
 *   RATIONAL ARRAYS   *
 **********************/

/*
 * Create and initialize an array of n rationals.
 */
rat_t* new_rat_array(const uint32_t n) {
  if (n > RAT_MAX_RATIONAL_ARRAY_SIZE) {
    EXIT("new_rational_array: requested array size is too large");
  }

  rat_t* a = (rat_t*)EGmalloc(n * sizeof(rat_t));
  rat_array_init(a, n);

  return a;
}

/*
 * Delete an array created by the previous function
 */
void free_rat_array(rat_t* a, const uint32_t n) {
  rat_array_clear(a, n);
  EGfree(a);
}
