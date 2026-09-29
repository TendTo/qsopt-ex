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
#ifndef RATIONALS_H
#define RATIONALS_H

#include <assert.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

/*
 * INTERNAL REPRESENTATION
 */

/*
 * A neorational is a union of size 64 bits.
 *
 * if the least bit is 1 it represents a
 * pointer to a gmp number.
 *
 * if the least bit is zero it is a struct consisting of a 32 bit
 * signed numerator, and a 31 bit unsigned denominator.
 */
typedef struct {
#ifdef WORDS_BIGENDIAN
  int32_t num;
  uint32_t den;
#else
  uint32_t den;
  int32_t num;
#endif
} qsx_rat32_t;

typedef struct {
  intptr_t gmp;
} qsx_ratgmp_t;

typedef union qsx_rational {
  qsx_rat32_t s;
  qsx_ratgmp_t p;
} rat_s;  // s for struct; p for pointer

typedef rat_s rat_t[1];
typedef rat_s* rat_ptr;

#define RAT_IS_RAT32 0x0
#define RAT_IS_RATGMP 0x1
/* the test for unicity in the denominator occurs frequently. the denominator is one if it is 2 :-) */
#define RAT_ONE_DEN 0x2

/*
 * Initialization: allocate and initialize
 * global variables.
 */
extern void init_rationals(void);

/*
 * Cleanup: free memory
 */
extern void cleanup_rationals(void);

/*
 * Set r to 0/1, Must be called before any operation on r.
 */
static inline void rat_init(rat_t r) {
  r->s.num = 0;
  r->s.den = RAT_ONE_DEN;
}

static inline void rat_init_pair(rat_t r1, rat_t r2) {
  rat_init(r1);
  rat_init(r2);
}

/*
 * Tests and conversions to/from gmp and rat32
 *
 * NOTE: the type mpq_ptr is defined in gmp.h. It's a pointer
 * to the internal structure representing a gmp number.
 */
static inline bool rat_is_rat32(const rat_t r) { return (r->p.gmp & RAT_IS_RATGMP) != RAT_IS_RATGMP; }

static inline bool rat_is_gmp(const rat_t r) { return (r->p.gmp & RAT_IS_RATGMP) == RAT_IS_RATGMP; }

static inline bool rat_is_uninitialized(const rat_t r) { return !rat_is_gmp(r) && r->s.den == 0; }

static inline mpq_ptr rat_get_mpq(const rat_t r) {
  assert(rat_is_gmp(r));
  return (mpq_ptr)(r->p.gmp ^ RAT_IS_RATGMP);
}

static inline mpz_ptr rat_get_mpq_num(const rat_t r) {
  assert(rat_is_gmp(r));
  return mpq_numref(rat_get_mpq(r));
}

static inline mpz_ptr rat_get_mpq_den(const rat_t r) {
  assert(rat_is_gmp(r));
  return mpq_denref(rat_get_mpq(r));
}

static inline int32_t rat_get_num(const rat_t r) {
  assert(rat_is_rat32(r));
  return r->s.num;
}

static inline uint32_t rat_get_den(const rat_t r) {
  assert(rat_is_rat32(r));
  return r->s.den >> 1;
}

static inline void rat_set_gmp(rat_t r, mpq_ptr gmp) { r->p.gmp = ((intptr_t)gmp) | RAT_IS_RATGMP; }

/*
 * Free mpq number attached to r if any, then set r to 0/1.
 * Must be called before r is deleted to prevent memory leaks.
 */
extern void rat_clear(rat_t r);

static inline void rat_clear_pair(rat_t r1, rat_t r2) {
  rat_clear(r1);
  rat_clear(r2);
}

/*
 * If r is represented as a gmp rational, convert it
 * to a pair of integers if possible.
 * If it's possible the gmp number is freed.
 */
extern void rat_normalize(rat_t r);

/*
 * ASSIGNMENT
 */

/*
 * Assign +1 or -1 to r
 */
static inline void rat_set_one(rat_t r) {
  rat_clear(r);
  r->s.num = 1;
}

static inline void rat_set_minus_one(rat_t r) {
  rat_clear(r);
  r->s.num = -1;
}

/*
 * Assignment operations: all set the value of the first argument (r or r1).
 * - in rat_set_int32 and rat_set_int64, num/den is normalized first
 *   (common factors are removed) and den must be non-zero
 * - in rat_set_mpq, q must be canonicalized first
 *   (a copy of q is made)
 * - rat_set copies r2 into r1 (if r2 is a gmp number,
 *   then a new gmp number is allocated with the same value
 *   and assigned to r1)
 * - rat_set_neg assigns the opposite of r2 to r1
 * - rat_set_abs assigns the absolute value of r2 to r1
 */
extern void rat_set_int32(rat_t r, int32_t num, uint32_t den);
extern void rat_set_int64(rat_t r, int64_t num, uint64_t den);
extern void rat_set32(rat_t r, int32_t num);
extern void rat_set64(rat_t r, int64_t num);
extern void rat_set_double(rat_t r, double d);

extern void rat_set_mpq(rat_t r, const mpq_t q);
extern void rat_set_mpq_prenormalize(rat_t r, const mpq_t q);
extern void rat_set_mpz(rat_t r, const mpz_t z);
extern void rat_set_mpf(rat_t r, const mpf_t f);
extern void rat_set(rat_t r1, const rat_t r2);
extern void rat_set_neg(rat_t r1, const rat_t r2);
extern void rat_set_abs(rat_t r1, const rat_t r2);

static inline void mpf_set_rat(mpf_t rop, const rat_t r) {
  if (rat_is_rat32(r)) {
    mpf_set_si(rop, rat_get_num(r));
    mpf_div_ui(rop, rop, rat_get_den(r));
  } else {
    mpf_set_q(rop, rat_get_mpq(r));
  }
}

static inline void mpq_set_rat(mpq_t rop, const rat_t r) {
  if (rat_is_rat32(r)) {
    mpq_set_si(rop, rat_get_num(r), rat_get_den(r));
  } else {
    mpq_set(rop, rat_get_mpq(r));
  }
}

/*
 * Copy r2 into r1: share the gmp index if r2
 * is a gmp number. Then clear r2.
 * This can be used without calling rat_init(r1).
 */
static inline void rat_copy_and_clear(rat_t r1, rat_t r2) {
  *r1 = *r2;
  r2->s.num = 0;
  r2->s.den = RAT_ONE_DEN;
}

/*
 * Swap values of r1 and r2
 */
static inline void rat_swap(rat_ptr r1, rat_ptr r2) {
  rat_s aux;

  aux = *r1;
  *r1 = *r2;
  *r2 = aux;
}

/*
 * Copy the numerator or denominator of r2 into r1
 */
extern void rat_set_num(rat_t r1, const rat_t r2);
extern void rat_set_den(rat_t r1, const rat_t r2);

/*
 * String parsing:
 * - set_from_string uses the GMP format with base 10:
 *       <optional_sign> <numerator>/<denominator>
 *       <optional_sign> <number>
 *
 * - set_from_string_base uses the GMP format with base b
 *
 * - set_from_float_string uses a floating point format:
 *   <optional sign> <integer part> . <fractional part>
 *   <optional sign> <integer part> <exp> <optional sign> <integer>
 *   <optional sign> <integer part> . <fractional part> <exp> <optional sign> <integer>
 *
 * where <optional sign> is + or - or nothing
 *       <exp> is either 'e' or 'E'
 *
 * The functions return -1 if the format is wrong and leave r unchanged.
 * The functions rat_set_from_string and rat_set_from_string_base return -2
 * if the denominator is 0.
 * Otherwise, the functions return 0 and the parsed value is stored in r.
 */
extern int rat_set_from_string(rat_t r, const char* s);
extern int rat_set_from_string_base(rat_t r, const char* s, int32_t b);
extern int rat_set_from_float_string(rat_t r, const char* s);

/*
 * ARITHMETIC: ALL OPERATIONS MODIFY THE FIRST ARGUMENT
 */

/*
 * Arithmetic operations:
 * - all operate on the first argument:
 *    rat_add: add r2 to r1
 *    rat_sub: subtract r2 from r1
 *    rat_mul: set r1 to r1 * r2
 *    rat_div: set r1 to r1/r2 (r2 must be nonzero)
 *    rat_neg: negate r
 *    rat_inv: invert r (r must be nonzero)
 *    rat_addmul: add r2 * r3 to r1
 *    rat_submul: subtract r2 * r3 from r1
 *    rat_addone: add 1 to r1
 *    rat_subone: subtract 1 from r1
 * - lcm/gcd operations (r1 and r2 must be integers)
 *    rat_lcm: store lcm(r1, r2) into r1
 *    rat_gcd: store gcd(r1, r2) into r1 (r1 and r2 must not be zero)
 * - floor and ceiling are also in-place operations:
 *    rat_floor: store largest integer <= r into r
 *    rat_ceil: store smaller integer >= r into r
 */
extern void rat_add(rat_t r1, const rat_t r2, const rat_t r3);
extern void rat_sub(rat_t r1, const rat_t r2, const rat_t r3);
extern void rat_mul(rat_t r1, const rat_t r2, const rat_t r3);
extern void rat_div(rat_t r1, const rat_t r2, const rat_t r3);
extern void rat_neg(rat_t r1, const rat_t r2);
extern void rat_inv(rat_t r1, const rat_t r2);
extern void rat_add_self(rat_t r1, const rat_t r2);
extern void rat_sub_self(rat_t r1, const rat_t r2);
extern void rat_mul_self(rat_t r1, const rat_t r2);
extern void rat_div_self(rat_t r1, const rat_t r2);
extern void rat_neg_self(rat_t r);
extern void rat_inv_self(rat_t r);

extern void rat_add_si(rat_t r1, const int64_t num, const uint64_t den);
extern void rat_sub_si(rat_t r1, const int64_t num, const uint64_t den);
extern void rat_mul_si(rat_t r1, const int64_t num, const uint64_t den);
extern void rat_div_si(rat_t r1, const int64_t num, const uint64_t den);

extern void rat_add_mpq(rat_t r1, const mpq_t q);
extern void rat_sub_mpq(rat_t r1, const mpq_t q);
extern void rat_mul_mpq(rat_t r1, const mpq_t q);
extern void rat_div_mpq(rat_t r1, const mpq_t q);

extern void rat_addmul(rat_t r1, const rat_t r2, const rat_t r3);
extern void rat_submul(rat_t r1, const rat_t r2, const rat_t r3);
extern void rat_add_one(rat_t r1);
extern void rat_sub_one(rat_t r1);

extern void rat_lcm(rat_t r1, const rat_t r2);
extern void rat_gcd(rat_t r1, const rat_t r2);

extern void rat_floor(rat_t r1, const rat_t r2);
extern void rat_ceil(rat_t r1, const rat_t r2);
extern void rat_floor_self(rat_t r);
extern void rat_ceil_self(rat_t r);

extern void rat_inv_mod(rat_t r, const rat_t mod);

/*
 * Exponentiation:
 * - rat_mulexp(r1, r2, n): multiply r1 by r2^n
 */
extern void rat_mulexp(rat_t r1, const rat_t r2, uint32_t n);

/*
 * Integer division and remainder
 * - r1 and r2 must both be integer
 * - r2 must be positive.
 * - Consider normalizing r2 before
 *
 * rat_integer_div(r1, r2) stores the quotient of r1 divided by r2 into r1
 * rat_integer_rem(r1, r2) stores the remainder into r1
 *
 * This implements the usual definition of division (unlike C).
 * If r = remainder and q = quotient then we have
 *    0 <= r < r2 and  r1 = q * r2 + r
 */
extern void rat_integer_div(rat_t r1, const rat_t r2);
extern void rat_integer_rem(rat_t r1, const rat_t r2);

/*
 * Generalized LCM: compute the smallest non-negative rational q
 * such that q/r1 is an integer and q/r2 is an integer.
 * - r1 and r2 can be arbitrary rationals.
 * - the result is stored in r1
 */
extern void rat_generalized_lcm(rat_t r1, const rat_t r2);

/*
 * Generalized GCD: compute the largest positive rational q
 * such that r1/q and r2/q are both integer.
 * - the result is stored in r1
 */
extern void rat_generalized_gcd(rat_t r1, const rat_t r2);

/*
 * SMT2 Versions of division and mod
 *
 * Intended semantics for div and mod:
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
extern void rat_smt2_div(rat_t q, const rat_t x, const rat_t y);
extern void rat_smt2_mod(rat_t q, const rat_t x, const rat_t y);

/*
 * TESTS: DO NOT MODIFY THE ARGUMENT(S)
 */

/*
 * Sign of r: rat_sgn(r) = 0 if r = 0
 *            rat_sgn(r) = +1 if r > 0
 *            rat_sgn(r) = -1 if r < 0
 */
static inline int rat_sgn(const rat_t r) {
  if (rat_is_gmp(r)) {
    return mpq_sgn(rat_get_mpq(r));
  } else {
    return (r->s.num < 0 ? -1 : (r->s.num > 0));
  }
}

/*
 * Compare r1 and r2:
 * - returns a negative number if r1 < r2
 * - returns 0 if r1 = r2
 * - returns a positive number if r1 > r2
 */
extern int rat_cmp(const rat_t r1, const rat_t r2);

/*
 * Compare r1 and num/den
 * - den must be nonzero
 * - returns a negative number if r1 < num/den
 * - returns 0 if r1 = num/den
 * - returns a positive number if r1 > num/den
 */
extern int rat_cmp_int32(const rat_t r1, int32_t num, uint32_t den);
extern int rat_cmp_int64(const rat_t r1, int64_t num, uint64_t den);

/*
 * Variants of rat_cmp:
 * - rat_le(r1, r2) is nonzero iff r1 <= r2
 * - rat_lt(r1, r2) is nonzero iff r1 < r2
 * - rat_gt(r1, r2) is nonzero iff r1 > r2
 * - rat_ge(r1, r2) is nonzero iff r1 >= r2
 * - rat_eq(r1, r2) is nonzero iff r1 = r2
 * - rat_neq(r1, r2) is nonzero iff r1 != r2
 */
static inline bool rat_le(const rat_t r1, const rat_t r2) { return rat_cmp(r1, r2) <= 0; }

static inline bool rat_lt(const rat_t r1, const rat_t r2) { return rat_cmp(r1, r2) < 0; }

static inline bool rat_ge(const rat_t r1, const rat_t r2) { return rat_cmp(r1, r2) >= 0; }

static inline bool rat_gt(const rat_t r1, const rat_t r2) { return rat_cmp(r1, r2) > 0; }

static inline bool rat_eq(const rat_t r1, const rat_t r2) { return rat_cmp(r1, r2) == 0; }

static inline bool rat_neq(const rat_t r1, const rat_t r2) { return rat_cmp(r1, r2) != 0; }

/*
 * Check whether r1 and r2 are opposite (i.e., r1 + r2 = 0)
 */
extern bool rat_opposite(const rat_t r1, const rat_t r2);

/*
 * Tests on rational r
 */
static inline bool rat_is_zero(const rat_t r) { return rat_is_gmp(r) ? mpq_sgn(rat_get_mpq(r)) == 0 : r->s.num == 0; }

static inline bool rat_is_nonzero(const rat_t r) {
  return rat_is_gmp(r) ? mpq_sgn(rat_get_mpq(r)) != 0 : r->s.num != 0;
}

static inline bool rat_is_one(const rat_t r) {
  return (r->s.den == RAT_ONE_DEN && r->s.num == 1) || (rat_is_gmp(r) && mpq_cmp_ui(rat_get_mpq(r), 1, 1) == 0);
}

static inline bool rat_is_minus_one(const rat_t r) {
  return (r->s.den == RAT_ONE_DEN && r->s.num == -1) || (rat_is_gmp(r) && mpq_cmp_si(rat_get_mpq(r), -1, 1) == 0);
}

static inline bool rat_is_pos(const rat_t r) { return (rat_is_rat32(r) ? r->s.num > 0 : mpq_sgn(rat_get_mpq(r)) > 0); }

static inline bool rat_is_nonneg(const rat_t r) {
  return (rat_is_rat32(r) ? r->s.num >= 0 : mpq_sgn(rat_get_mpq(r)) >= 0);
}

static inline bool rat_is_neg(const rat_t r) { return (rat_is_rat32(r) ? r->s.num < 0 : mpq_sgn(rat_get_mpq(r)) < 0); }

static inline bool rat_is_nonpos(const rat_t r) {
  return (rat_is_rat32(r) ? r->s.num <= 0 : mpq_sgn(rat_get_mpq(r)) <= 0);
}

static inline bool rat_is_integer(const rat_t r) {
  return (rat_is_rat32(r) && r->s.den == RAT_ONE_DEN) || (rat_is_gmp(r) && mpz_cmp_ui(mpq_denref(rat_get_mpq(r)), 1UL) == 0);
}

/*
 * DIVISIBILITY TESTS
 */

/*
 * Check whether r1 divides r2: both must be integers
 * - r1 must be non-zero
 * - side effect: r1 is normalized
 */
extern bool rat_integer_divides(rat_t r1, const rat_t r2);

/*
 * General divisibility check: return true iff r2/r1 is an integer
 * - r1 must be non-zero
 */
extern bool rat_divides(const rat_t r1, const rat_t r2);

/*
 * SMT2 version:
 * - if t1 is non-zero, return true iff (r2/r1) is an integer
 * - if t1 is zero, return true iff r2 is zero
 */
extern bool rat_smt2_divides(const rat_t r1, const rat_t r2);

/*
 * Tests on integer rational r mod m
 */
static inline bool rat_is_zero_mod(const rat_t r, const rat_t m) {
  assert(rat_is_integer(r) && rat_is_integer(m) && rat_is_pos(m));
  return rat_integer_divides((rat_ptr)m, r);
}

static inline bool rat_is_nonzero_mod(const rat_t r, const rat_t m) { return !rat_is_zero_mod(r, m); }

/*
 * CONVERSIONS TO INTEGERS
 */

/*
 * Check whether r is a small integer (use with care: this
 * ignores the case where r is a small integer that happens
 * to be represented as a gmp rational).
 * Call rat_normalize(r) first if there's a doubt.
 */
static inline bool rat_is_smallint(rat_t r) { return r->s.den == RAT_ONE_DEN; }

/*
 * Convert r to an integer, provided rat_is_smallint(r) is true
 */
static inline int32_t rat_get_smallint(rat_t r) { return rat_get_num(r); }

/*
 * Conversions: all functions attempt to convert r into an integer or
 * a pair of integers (num/den). If the conversion is not possible
 * the functions return false. Otherwise, the result is true and the
 * value is returned in v or num/den.
 */
extern bool rat_get32(rat_t r, int32_t* v);
extern bool rat_get64(rat_t r, int64_t* v);
extern bool rat_get_int32(rat_t r, int32_t* num, uint32_t* den);
extern bool rat_get_int64(rat_t r, int64_t* num, uint64_t* den);

/*
 * Similar to the conversion functions above but just check
 * whether the conversions are possible.
 */
extern bool rat_is_int32(const rat_t r);    // r is a/1 where a is int32
extern bool rat_is_int64(const rat_t r);    // r is a/1 where a is int64
extern bool rat_fits_int32(const rat_t r);  // r is a/b where a is int32, b is uint32
extern bool rat_fits_int64(const rat_t r);  // r is a/b where a is int64, b is uint64

/*
 * Size estimate
 * - this returns approximately the number of bits to represent r's numerator
 * - this may not be exact (typically rounded up to a multiple of 32)
 * - also if r is really really big, this function may return
 *   UINT32_MAX (not very likely!)
 */
extern uint32_t rat_size(const rat_t r);

/*
 * CONVERSION TO GMP OBJECTS
 */

/*
 * Store r into the GMP integer z.
 * - return false if r is not a integer, true otherwise
 */
extern bool rat_to_mpz(const rat_t r, mpz_t z);

extern void rat_to_mpq(const rat_t r, mpq_t q);

/*
 * Convert to a floating point number
 */
extern double rat_to_d(const rat_t r);

/*
 * Set from floating point numbers
 */
static inline void rat_set_float(rat_t r, float val) { rat_set_double(r, (double)val); }

/*
 * PRINT
 */

/*
 * Print r on stream f.
 * rat_print_abs prints the absolute value
 */
extern void rat_print(FILE* f, const rat_t r);
extern void rat_print_abs(FILE* f, const rat_t r);

/*
 * Return a string representation of r.
 * The string is allocated on the heap and must be freed by the caller.
 */
extern char* rat_get_str(const rat_t r);

/*
 * HASH FUNCTION
 */

/*
 * Hash functions: return a hash of numerator or denominator.
 * - hash_numerator(r) = numerator MOD bigprime
 * - hash_denominator(r) = denominator MOD bigprime
 * where bigprime is the largest prime number smaller than 2^32.
 */
extern uint32_t rat_hash_numerator(const rat_t r);
extern uint32_t rat_hash_denominator(const rat_t r);
extern void rat_hash_decompose(const rat_t r, uint32_t* h_num, uint32_t* h_den);

/*
 * RATIONAL ARRAYS
 */

#define RAT_MAX_RATIONAL_ARRAY_SIZE (UINT32_MAX / sizeof(rat_t))

/*
 * Create and initialize an array of n rationals.
 */
extern rat_t* new_rat_array(uint32_t n);

/*
 * Delete an array created by the previous function
 * - n must be the size of a
 */
extern void free_rat_array(rat_t* a, uint32_t n);

/*
 * Initialize an array of n rationals.
 * - The array must have been allocated first.
 * - n must be the size of a
 */
static inline void rat_array_init(rat_t* a, uint32_t n) {
  for (uint32_t i = 0; i < n; i++) {
    rat_init(a[i]);
  }
}

/*
 * Clear an array of n rationals.
 * The array itself is not freed.
 * - n must be the size of a
 */
static inline void rat_array_clear(rat_t* a, uint32_t n) {
  for (uint32_t i = 0; i < n; i++) {
    rat_clear(a[i]);
  }
}

#endif
