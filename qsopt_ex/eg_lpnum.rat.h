/* EGlib "Efficient General Library" provides some basic structures and
 * algorithms commons in many optimization algorithms.
 *
 * Copyright (C) 2005 Daniel Espinoza and Marcos Goycoolea.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA 
 * */
#ifndef __EG_LPNUM_RAT__
#define __EG_LPNUM_RAT__

#include <stdlib.h>
#include <string.h>

#include <gmp.h>
#include "rationals.h"

#include "eg_lpnum.h"

/** @file
 * @brief Interface for rational implementation of the EGlpNum_t type.
 * @par History
 * - 2006-02-01
 * 					- Add verbosity flag for continued fraction conversions.
 * @ingroup EGlpNum */
/** @addtogroup EGlpNum */
/** @{ */
/* ========================================================================= */
/** extern definitions of constaants for different set-ups */
extern const rat_t __zeroLpNum_rat__;
extern const rat_t __oneLpNum_rat__;
extern const rat_t __MaxLpNum_rat__;
extern const rat_t __MinLpNum_rat__;
#define rat_zeroLpNum __zeroLpNum_rat__
#define rat_oneLpNum  __oneLpNum_rat__
#define rat_epsLpNum  __zeroLpNum_rat__
#define rat_MaxLpNum  __MaxLpNum_rat__
#define rat_MinLpNum  __MinLpNum_rat__

/* ========================================================================= */
/** @brief This function read a number in float form and store it in an rat_t
 * variable, returning how many chars read to create the number, the twist is
 * that it does an 'exact' transformation, in the sense that 0.33333333 will be
 * stored as 33333333/100000000.
 * @param str input string.
 * @param var variable where we will store the number as rational.
 * @return number of reade chars.
 * @par Descriptiom:
 * If the input string doesn't contain a number, 'var' will be set to zero, and
 * the number of readed chars will be zero, we assume that there are no leading
 * empty spaces, (nor tabs), no eschape characters, and trailing zeros are
 * considered as readed, but no the following spaces. 
 * @note
 * This function will only read number in decimal base, in a future release we
 * may include a more general reader. */
int rat_EGlpNumReadStrXc (rat_t var,
													char const *str);

/* ========================================================================= */
/** @brief Read from a string a number and store it in the given rat_t, return
 * the number of chars readed from the input string */
#define rat_EGlpNumReadStr(a,str) rat_EGlpNumReadStrXc(a,str)

/* ========================================================================= */
/** @brief given a rat_t, write it to a string (to be allocated internally), 
 * and return it. */
#define rat_EGlpNumGetStr(a) rat_get_str(a);

/* ========================================================================= */
/** @brief given an array of type rat_t, free it, if the pointer is NULL
 * nothing happen. */
#define rat_EGlpNumFreeArray(ea) ({\
	size_t __sz = __EGlpNumArraySize(ea);\
	rat_t* __ptr__ = (ea);\
	while(__sz--) rat_clear(__ptr__[__sz]);\
	__EGlpNumFreeArray(ea);})


/* ========================================================================= */
/** @brief Reallocate and initialize (if needed) 'size' elements of type 
 * rat_t and return it, if no more memory, exit(1) */
#define rat_EGlpNumReallocArray(lptr, lsize) ({\
	rat_t **__ptr__ = (lptr);\
	size_t *__ntmp__ = (size_t *) *__ptr__, __sz__ = (lsize);\
	size_t __psz__;\
	/* if no memory allocated before we just call the regular allocator */\
	if (!*__ptr__)\
		*__ptr__ = rat_EGlpNumAllocArray (__sz__);\
	else\
	{\
		/* first check that the previous size is not larger than the current */\
		__ntmp__--;\
		__psz__ = __ntmp__[0];\
		if (__psz__ < __sz__)\
		{\
			/* now we have to do the reallocation */\
			*__ptr__ = (rat_t *) __ntmp__;\
			*__ptr__ = EGrealloc(*__ptr__,sizeof(rat_t) * __sz__ + sizeof(size_t));\
			__ntmp__ = (size_t *) *__ptr__;\
			__ntmp__[0] = __sz__;\
			__ntmp__++;\
			*__ptr__ = (rat_t *) __ntmp__;\
			for (; __psz__ < __sz__; __psz__++) rat_init ((*__ptr__)[__psz__]);\
		}\
	}\
})

/* ========================================================================= */
/** @brief Allocate and initialize (if needed) 'size' elements of type rat_t
 * and return it, if no more memory, exit(1) */
#define rat_EGlpNumAllocArray(size) ({\
	size_t __i__ = (size);\
	rat_t *__res = __EGlpNumAllocArray(rat_t,__i__);\
	while(__i__--) rat_init(__res[__i__]);\
	__res;})

/* ========================================================================= */
/** @brief set the given rational number , to the value to the value of the 
 * given mpf_t, this conversion is done using the continuous fraction method.
 * @param var rat_t where we will store the value.
 * @param flt mpf_t value to be stored in 'var'.
 * @par Description:
 * This function is intended to set initial values to variables. Note also 
 * that if the number is
 * writen in the form \f$x=\bar{x}\cdot 2^e\f$ with \f$0.5<|\bar{x}|<1\f$, 
 * then \f$\left|x-\frac{p}{q}\right|<2^{e-EGLPNUM_PRECISION}\f$.
 * */
void rat_EGlpNumSet_mpf (rat_t var,
												 mpf_t flt);

/* ========================================================================= */
/** @brief set the given number pointer, set its value to the given double.
 * @param var rat_t where we will store the double value.
 * @param dbl double value to be stored in 'var'.
 * @par Description:
 * This function is intended to set initial values to variables; note that the
 * double is a number and not a pointer to that value, be carefull with this
 * detail. Also, due to implementation details this function can't deal with
 * numbers above 1e158 or smaller than 1e-158. Note also that if the number is
 * writen in the form \f$x=\bar{x}\cdot 2^e\f$ with \f$0.5<|\bar{x}|<1\f$, 
 * then \f$\left|x-\frac{p}{q}\right|<2^{e-64}\f$.
 * */
void rat_EGlpNumSet (rat_t var,
										 const double dbl);

/* ========================================================================= */
/** @brief Stores in the first number the ceil value of the second number, i.e.
 * EGlpNumCeil(a,b) <==> a= ceil(b) */
#define rat_EGlpNumCeil(a, b) rat_ceil(a, b)

/* ========================================================================= */
/** @brief Stores in the first number the floor value of the second number, i.e.
 * EGlpNumFloor(a,b) <==> a= floor(b) */
#define rat_EGlpNumFloor(a, b) rat_floor(a, b)

/* ========================================================================= */
/** @brief store the (multiplicative) inverse of a number to itself, i.e.
 * implement a = 1/a.
 * @param a the number to be inverted. */
#define rat_EGlpNumInv(a) rat_inv(a,a)

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a rat_t first number to compare.
 * @param b rat_t second number to compare.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' return 1 if a == b, otherwise it return 0
 * */
#define rat_EGlpNumIsEqqual(a,b) (rat_eq(a,b))

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a rat_t first number to compare.
 * @param b rat_t second number to compare.
 * @param error rat_t maximum difference allowed between both
 * numbers.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' and a tolerance 'error',
 * return 1 if |a-b|<= error, otherwise it return 0.
 * */
#define rat_EGlpNumIsEqual(a,b,error) (rat_eq(a,b))
#define rat_EGlpNumIsNeq(a,b,error)   (!(rat_eq(a,b)))
#define rat_EGlpNumIsNeqq(a,b)        (!(rat_eq(a,b)))
#define rat_EGlpNumIsNeqZero(a,error) (rat_sgn(a))
#define rat_EGlpNumIsNeqqZero(a)      (rat_sgn(a))

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a rat_t the first number.
 * @param b rat_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define rat_EGlpNumIsLess(a,b) (rat_cmp(a,b) < 0)

/* ========================================================================= */
/** @brief test if a number is greater than zero
 * @param a number to test
 * @return int one if success, zero otherwise.
 * */
#define rat_EGlpNumIsGreatZero(a) (rat_sgn(a) > 0)

/* ========================================================================= */
/** @brief test if a number is less than zero
 * @param a number to test
 * @return int one if success, zero otherwise.
 * */
#define rat_EGlpNumIsLessZero(a) (rat_sgn(a) < 0)

/* ========================================================================= */
/** @brief test if the sum of the first two numbers is less thatn the third
 * number.
 * @param a rat_t the first number.
 * @param b rat_t the second number
 * @param c rat_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a + b < c), zero toherwise.
 * */
#define rat_EGlpNumIsSumLess(a, b, c) ({\
	rat_t __lpnum__;int __res__=0;rat_init(__lpnum__);\
	rat_add(__lpnum__,a,b);\
	__res__=(rat_cmp (__lpnum__, c) < 0);\
	rat_clear(__lpnum__);\
	__res__;\
})

/* ========================================================================= */
/** @brief test if the diference of the first two numbers is less thatn the 
 * third number.
 * @param a rat_t the first number.
 * @param b rat_t the second number
 * @param c rat_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a - b < c), zero toherwise.
 * */
#define rat_EGlpNumIsDiffLess(a, b, c) ({\
	rat_t __lpnum__;int __res__=0;rat_init(__lpnum__);\
	rat_sub (__lpnum__, a, b);\
	__res__=(rat_cmp (__lpnum__, c) < 0);\
	rat_clear(__lpnum__);\
	__res__;\
})

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a rat_t the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define rat_EGlpNumIsLessDbl(a,b) (rat_to_d(a) < b)

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a rat_t the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a > b, zero
 * otherwise.
 * */
#define rat_EGlpNumIsGreaDbl(a,b) (rat_to_d(a) > b)

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a rat_t the first number.
 * @param b rat_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a <= b, zero
 * otherwise.
 * */
#define rat_EGlpNumIsLeq(a,b) (rat_cmp(a,b) <= 0)

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a rat_t source number (it won't change value).
 * @param b rat_t source number (it won't change value).
 * @param c rat_t denominator of the difference (it won't change value).
 * @param d rat_t where to store the value .
 * @par Description:
 * Set @f$a = \frac{b - c}{d} @f$ */
#define rat_EGlpNumCopyDiffRatio(a, b, c, d) ({\
	rat_sub (a, b, c);\
	rat_div (a, a, d);\
})

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a rat_t source number (it won't change value).
 * @param b rat_t source number (it won't change value).
 * @param dest rat_t where to store the value stored in 'orig'.
 * @par Description:
 * Set dest = a - b */
#define rat_EGlpNumCopyDiff(dest,a,b) rat_sub(dest,a,b)

/* ========================================================================= */
/** @brief copy the value of the sum of the second and third parameter
 * @param a rat_t source number (it won't change value).
 * @param b rat_t source number (it won't change value).
 * @param dest rat_t where to store the sum.
 * @par Description:
 * Set dest = a + b */
#define rat_EGlpNumCopySum(dest,a,b) rat_add(dest,a,b)

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param orig rat_t source number (it won't change value).
 * @param dest rat_t where to store the value stored in 'orig'.
 * @par Description:
 * Given two numbers copy the values in 'orig', into 'dest'.
 * */
#define rat_EGlpNumCopy(dest,orig) rat_set(dest,orig)

/* ========================================================================= */
/** @brief change the fist number to the maximum between itself and the 
 * absolute value of the second.
 * @param orig rat_t source number (it won't change value).
 * @param dest rat_t where to store the value stored in 'orig'.
 * @par Description:
 * implement dest = max(dest,abs(orig))
 * */
#define rat_EGlpNumSetToMaxAbs(dest, orig) ({\
	rat_t __lpnum__;rat_init(__lpnum__);\
	rat_set_abs (__lpnum__, orig);\
	if (rat_cmp (dest, __lpnum__) < 0)\
		rat_set (dest, __lpnum__);\
	rat_clear(__lpnum__);\
})

#define rat_EGlpNumSetToMinAbs(dest, orig) ({\
	rat_t __lpnum__;rat_init(__lpnum__);\
	rat_set_abs (__lpnum__, orig);\
	if (rat_cmp (dest, __lpnum__) > 0)\
		rat_set (dest, __lpnum__);\
	rat_clear(__lpnum__);\
})

/* ========================================================================= */
/** @brief copy the square of the second argument, divided by the third 
 * argument into the first argument.
 * @param dest rat_t where to store the result
 * @param orig rat_t second parameter
 * @param den rat_t third parameter
 * @par Description:
 * compute dest = (orig*orig)/den
 * */
#define rat_EGlpNumCopySqrOver(dest, orig, den) ({\
	rat_mul (dest, orig, orig);\
	rat_div (dest, dest, den);\
})

/* ========================================================================= */
/** @brief copy the value of the absolute value of the second parameter to the 
 * first parameter.
 * @param orig rat_t source number (it won't change value).
 * @param dest rat_t where to store the absolute value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy its absolute value to 'dest'. i.e.
 * dest = |orig|
 * */
#define rat_EGlpNumCopyAbs(dest,orig) rat_set_abs(dest,orig)

/* ========================================================================= */
/** @brief copy minus the value of the second parameter to the 
 * first parameter.
 * @param orig rat_t the source number (it won't change value).
 * @param dest rat_t where to store minus the value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy minus the value to 'dest'. i.e.
 * dest = -orig
 * */
#define rat_EGlpNumCopyNeg(dest,orig) rat_neg(dest,orig)

/* ========================================================================= */
/** @brief Set des = op1/op2.
 * @param dest rat_t where we will store the result.
 * @param op1 rat_t numerator of the fraction (possibly non an integer)
 * @param op2 rat_t denominator of the fraction (possibly non an integer)
 * @par Description:
 *  Set des = op1/op2
 * */
#define rat_EGlpNumCopyFrac(dest,op1,op2) rat_div(dest,op1,op2)

/* ========================================================================= */
/** @brief copy the first 'size' values in the second array to the first array.
 * @param orig rat_t* pointer to the array from where we will copy the
 * values (it won't change value).
 * @param dest rat_t* pointer to where to store the first 'size' values 
 * stored in 'orig'.
 * @param size unsigned int specifying how many values of 'orig' will be copied
 * onto 'dest'
 * @par Description:
 * This function is provided to (possible) make fast copies of arrays of
 * numbers, the arrays should be of length at least 'size', and the resulting
 * copy is absolutely independent froom the original, any change in one vale of
 * one array won't change values on the other array.
 * */
#define rat_EGlpNumCopyArray(dest,orig,size) {\
	register unsigned int __i__ = size;\
	for(;__i__--;)\
	{\
		rat_set(dest[__i__],orig[__i__]);\
	}\
}

/* ========================================================================= */
/** @brief Sub to a given number the product of two numbers.
 * @param a rat_t the number that we are going to Sub to.
 * @param b rat_t value to be multiplyed.
 * @param c rat_t value to be multiplyed.
 * @par Description:
 * This function implements a = a - b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define rat_EGlpNumSubInnProdTo(a, b, c) ({\
	rat_t __lpnum__;rat_init(__lpnum__);\
	rat_mul (__lpnum__, b, c);\
	rat_sub (a, a, __lpnum__);\
	rat_clear(__lpnum__);\
})

/* ========================================================================= */
/** @brief Add to a given number the product of two numbers.
 * @param a rat_t the number that we are going to add to.
 * @param b rat_t value to be multiplyed.
 * @param c rat_t value to be multiplyed.
 * @par Description:
 * This function implements a = a + b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define rat_EGlpNumAddInnProdTo(a, b, c) ({\
	rat_t __lpnum__;rat_init(__lpnum__);\
	rat_mul (__lpnum__, b, c);\
	rat_add (a, a, __lpnum__);\
	rat_clear(__lpnum__);\
})

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a rat_t the number that we are going to substract to.
 * @param b unsigned int value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define rat_EGlpNumSubUiTo(a,b) ({ \
	rat_t __lpnum__;rat_init(__lpnum__);\
	rat_set64(__lpnum__,(long int)b);\
	rat_sub_self(a, __lpnum__); \
	rat_clear(__lpnum__); \
})

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a rat_t the number that we are going to add to.
 * @param b unsigned int value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define rat_EGlpNumAddUiTo(a,b) ({ \
	rat_t __lpnum__;rat_init(__lpnum__);\
	rat_set64(__lpnum__,(long int)b);\
	rat_add_self(a, __lpnum__); \
	rat_clear(__lpnum__); \
})

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a rat_t the number that we are going to add to.
 * @param b rat_t value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define rat_EGlpNumAddTo(a,b) rat_add(a,a,b)

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a rat_t the number that we are going to substract
 * from.
 * @param b rat_t value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define rat_EGlpNumSubTo(a,b) rat_sub(a,a,b)

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a rat_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b rat_t value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define rat_EGlpNumMultTo(a,b) rat_mul(a,a,b)

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a rat_t the number that we are going to divide by
 * the second number and store the result.
 * @param b rat_t value to be divide to 'a'.
 * @par Description:
 * This function implements a = a / b, and clearly don't change the value
 * stored in 'b'.
 * */
#define rat_EGlpNumDivTo(a,b) rat_div(a,a,b)

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a rat_t the number that we are going to divide by
 * the second number and store the result.
 * @param b unsigned int value to be divided to 'a'.
 * @par Description:
 * This function implements a = a / b, and don't change the value
 * stored in 'b'.
 * */
#define rat_EGlpNumDivUiTo(a,b) do{ \
	rat_t __lpnum__;rat_init(__lpnum__);\
	rat_set64(__lpnum__,(long int)b);\
	rat_div_self(a, __lpnum__); \
	rat_clear(__lpnum__); \
}while(0)

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a rat_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b unsigned int value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define rat_EGlpNumMultUiTo(a,b) do{ \
	rat_t __lpnum__;rat_init(__lpnum__);\
	rat_set64(__lpnum__,(long int)b);\
	rat_mul_self(a, __lpnum__); \
	rat_clear(__lpnum__); \
}while(0)

/* ========================================================================= */
/** @brief Reset the value of the pointed number to zero.
 * @param a rat_t the value to be set to zero.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 0;
 * */
#define rat_EGlpNumZero(a) rat_clear(a)

/* ========================================================================= */
/** @brief Reset the value of the pointed number to one.
 * @param a rat_t value to be set to one.
 * @par Descrpition:
 * Reset a to one, i.e. implements a = 1;
 * */
#define rat_EGlpNumOne(a) rat_set_one(a)

/* ========================================================================= */
/** @brief Change the sign of the number.
 * @param a rat_t number we will change sign.
 * @par Descrpition:
 * Change the sign of the given number, i.e. implements a = -a
 * */
#define rat_EGlpNumSign(a) rat_neg(a,a)

/* ========================================================================= */
/** @brief return the closest double value of the given pointer number.
 * @param a rat_t number that we will be transformed to double.
 * @return double the closest double representation of the given number.
 * par Description:
 * return the double number closest in value to the value stored in a.
 * */
#define rat_EGlpNumToLf(a) rat_to_d(a)

/* ========================================================================= */
/** @brief initialize the internal memory of a given variable */
#define rat_EGlpNumInitVar(a) rat_init(a)

/* ========================================================================= */
/** @brief free the internal memory of a given variable */
#define rat_EGlpNumClearVar(a) rat_clear(a)

/* ========================================================================= */
/** @} */
#endif
