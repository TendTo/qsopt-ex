
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "qs_config.h"
#include "logging-private.h"

#include "QSopt_ex.h"
/* ========================================================================= */
/** @name Static Variables
 * Set of static variables for this program. */
/*@{*/
/** @brief range of values for first OA level */
static unsigned s1 = 3;
/** @brief range of values for second OA level */
static unsigned s2 = 3;
/** @brief number of columns for first OA level */
static unsigned k1 = 3;
/** @brief number of columns for second OA level */
static unsigned k2 = 3;
/** @brief use (or not) scaling of the resulting LP */
static int use_scaling = 1;
/** @brief use dlouble floating point output */
static int use_double = 0;
/** @brief output feasibility LP */
static int use_feas = 0;
/** @brief use MPS file format */
static int use_mps = 0;
/** @brief file name output. */
const char *out_file = NULL;
/** @brief Strength level required */
static unsigned t = 2;
/** @brief temporal string space */
char strtmp[1024];
/** @brief whether to output the dual problem */
static int use_dual = 0;
/*@}*/

/* ========================================================================= */
/** @brief Display options to the screen */
static inline void sl_usage (char const *const s)
{
	fprintf (stderr,
					 "This programs compute bounds for mixed-level orthogonal arrays (OA) of different strengths for two levels\n");
	fprintf (stderr, "Usage: %s [options]\n", s);
	fprintf (stderr, "    -a n   range of values for the first OA level (>=2)\n");
	fprintf (stderr,
					 "    -b n   range of values for the second OA level (>=2)\n");
	fprintf (stderr,
					 "    -c n   number of columns for the first OA level (>=2)\n");
	fprintf (stderr,
					 "    -e n   number of columns for the second OA level (>=2)\n");
	fprintf (stderr,
					 "    -f n   if n > 0, output the feasibility problem (zero objective), otherwise, the optimality problem\n");
	fprintf (stderr,
					 "    -d n   if n > 0, write the LP in double precition arithmetic, otherwise, write it in rational form\n");
	fprintf (stderr, "    -o f   filename where to write the output LP\n");
	fprintf (stderr,
					 "    -s n   if n > 0 scale the resulting LP, otherwise present the unscaled LP\n");
	fprintf (stderr, "    -t n   required strength of the OA (>=1)\n");
	fprintf (stderr, "    -m n   if n > 0 output in MPS format, otherwise, write it in LP format\n");
	fprintf (stderr, "    -i n   if n > 0 output the dual problem, otherwise, output the primal problem\n");
}

/* ========================================================================= */
/** @brief Display options to the screen */
static inline int sl_parseargs (int argc,
																char **argv)
{
	int c;
	while ((c = getopt (argc, argv, "a:b:c:e:f:m:s:d:o:t:i:")) != EOF)
	{
		switch (c)
		{
		case 'a':
			s1 = atoi (optarg);
			break;
		case 'b':
			s2 = atoi (optarg);
			break;
		case 'c':
			k1 = atoi (optarg);
			break;
		case 'e':
			k2 = atoi (optarg);
			break;
		case 's':
			use_scaling = atoi (optarg);
			break;
		case 'd':
			use_double = atoi (optarg);
			break;
		case 'f':
			use_feas = atoi (optarg);
			break;
		case 'm':
			use_mps = atoi (optarg);
			break;
		case 'o':
			out_file = optarg;
			break;
		case 't':
			t = atoi (optarg);
			break;
		case 'i':
			use_dual = atoi (optarg);
			break;
		default:
			printf ("Unknown option %c\n", c);
			sl_usage (argv[0]);
			return 1;
		}
	}
	if (s1 < 2 || s2 < 2 || k1 < 1 || k2 < 1 || t < 1)
	{
		sl_usage (argv[0]);
		return 1;
	}
	if (!out_file)
	{
		if (use_mps)
			out_file = "output.mps";
		else
			out_file = "output.lp";
	}
	fprintf (stderr, "Running %s\nOptions:\n", argv[0]);
	fprintf (stderr, "\tfirst level columns = %d\n", k1);
	fprintf (stderr, "\tsecond level columns = %d\n", k2);
	fprintf (stderr, "\tfirst level range = %d\n", s1);
	fprintf (stderr, "\tsecond level range = %d\n", s2);
	fprintf (stderr, "\tt = %d\n", t);
	fprintf (stderr, "\toutputting %s problem\n", use_feas ? "feasibility" : "optimality");
	fprintf (stderr, "\t%s scaling\n", use_scaling ? "using" : "not using");
	fprintf (stderr, "\t%s output\n", use_double ? "double" : "rational");
	fprintf (stderr, "\toutput file : %s\n", out_file);
	return 0;
}

/* ========================================================================= */
/** @brief compute the krawtchouk ppolynomial, defined as
 * \f[P^s_j(x,m)=\sum\limits_{i=0}^j(-1)^i(s-1)^{j-i}\binom{x}{i}\binom{m-x}{j-i}\f]
 * @param x the \f$x\f$ parameter of the function.
 * @param s the \f$s\f$ parameter of the function.
 * @param j the \f$j\f$ parameter of the function.
 * @param m the \f$m\f$ parameter of the function.
 * @param rop where to store the result */
void Kpoly (unsigned x,
						unsigned s,
						unsigned j,
						unsigned m,
						mpz_t rop)
{
	register unsigned i = j + 1;
	mpz_t z1,
	  z2;
	mpz_init (z1);
	mpz_init (z2);
	mpz_set_ui (rop, (unsigned long)0);
	EXIT (m < x, "m < x! impossible!");
	EXIT (s < 1, "s < 1! impossible!");
	while (i--)
	{
		mpz_bin_uiui (z1, (unsigned long)x, (unsigned long)i);
		mpz_set (z2, z1);
		mpz_bin_uiui (z1, (unsigned long)(m - x), (unsigned long)(j - i));
		mpz_mul (z2, z2, z1);
		mpz_ui_pow_ui (z1, (unsigned long)(s - 1), (unsigned long)(j - i));
		mpz_mul (z2, z2, z1);
		if (i & 1U)
			mpz_sub (rop, rop, z2);
		else
			mpz_add (rop, rop, z2);
	}
	mpz_clear (z1);
	mpz_clear (z2);
}

mpq_QSdata * to_dual(mpq_QSdata *p_mpq, char problem_name[]) {
	mpq_QSdata * dual_p_mpq = 0;
	mpq_t oneLpNum;
	mpq_init(oneLpNum);
	mpq_set_si(oneLpNum, 1, 1);

	dual_p_mpq = mpq_QScreate_prob (problem_name, QS_MIN);
	int n_cols = mpq_QSget_colcount(p_mpq);
	int n_rows = mpq_QSget_rowcount(p_mpq);

	mpq_t obj_coeffs[n_cols];
	for (int i = 0; i < n_cols; i++) {
		mpq_init(obj_coeffs[i]);
	}
	mpq_QSget_obj(p_mpq, obj_coeffs);
	mpq_t rhss[n_rows];
	for (int i = 0; i < n_rows; i++) {
		mpq_init(rhss[i]);
	}
	mpq_QSget_rhs(p_mpq, rhss);

	char senses[n_rows];
	mpq_QSget_senses(p_mpq, senses);

	for (int row = 0; row < n_rows; row++) {
		if (senses[row] == 'E') {
			mpq_QSnew_col(dual_p_mpq, rhss[row], mpq_NINFTY, mpq_INFTY, NULL);
		} else if (senses[row] == 'G') {
			mpq_QSnew_col(dual_p_mpq, rhss[row], mpq_NINFTY, mpq_zeroLpNum, NULL);
		} else if (senses[row] == 'L') {
			mpq_QSnew_col(dual_p_mpq, rhss[row], mpq_zeroLpNum, mpq_INFTY, NULL);
		} else {
			fprintf(stderr, "Unknown sense %c for row %d\n", senses[row], row);
			exit(1);
		}
	}

	for (int col = 0; col < n_cols; col++) {
		mpq_t lower, upper;
		mpq_init(lower);
		mpq_init(upper);

		mpq_QSget_bound(p_mpq, col, 'L', &lower);
		mpq_QSget_bound(p_mpq, col, 'U', &upper);

		if (mpq_cmp(lower, mpq_zeroLpNum) != 0 && mpq_cmp(lower, mpq_NINFTY) != 0) {
			mpq_QSnew_col(dual_p_mpq, lower, mpq_NINFTY, mpq_zeroLpNum, NULL);
			printf("Column %d has lower bound %s\n", col, mpq_get_str(NULL, 10, lower));
		} else if (mpq_cmp(upper, mpq_zeroLpNum) != 0 && mpq_cmp(upper, mpq_INFTY) != 0) {
			mpq_QSnew_col(dual_p_mpq, upper, mpq_zeroLpNum, mpq_INFTY, NULL);
			printf("Column %d has upper bound %s\n", col, mpq_get_str(NULL, 10, upper));
		} else {
			continue;
		}

	}

	for (int col = 0; col < n_cols; col++) {
		mpq_t lower, upper;
		mpq_init(lower);
		mpq_init(upper);

		mpq_QSget_bound(p_mpq, col, 'L', &lower);
		mpq_QSget_bound(p_mpq, col, 'U', &upper);

		char sense = '?';
		if (mpq_cmp(lower, mpq_zeroLpNum) == 0) {
			mpq_QSnew_row(dual_p_mpq, obj_coeffs[col], 'L', NULL);
		} else if (mpq_cmp(upper, mpq_zeroLpNum) == 0) {
			mpq_QSnew_row(dual_p_mpq, obj_coeffs[col], 'G', NULL);
		} else if (mpq_cmp(lower, mpq_NINFTY) == 0 && mpq_cmp(upper, mpq_INFTY) == 0) {
			mpq_QSnew_row(dual_p_mpq, obj_coeffs[col], 'E', NULL);
		} else if (mpq_cmp(lower, mpq_zeroLpNum) != 0 && mpq_cmp(lower, mpq_NINFTY) != 0) {
			mpq_QSnew_row(dual_p_mpq, obj_coeffs[col], 'E', NULL);
			mpq_QSchange_coef(dual_p_mpq, col, mpq_QSget_colcount(dual_p_mpq) - 1, oneLpNum);
		} else if (mpq_cmp(upper, mpq_zeroLpNum) != 0 && mpq_cmp(upper, mpq_INFTY) != 0) {
			mpq_QSnew_row(dual_p_mpq, obj_coeffs[col], 'E', NULL);
			mpq_QSchange_coef(dual_p_mpq, col, mpq_QSget_colcount(dual_p_mpq) - 1, oneLpNum);
		} else {
			printf("Column %d has lower bound %s and upper bound %s\n", col, mpq_get_str(NULL, 10, lower), mpq_get_str(NULL, 10, upper));
			exit(1);
		}

		mpq_clear(lower);
		mpq_clear(upper);

		for (int row = 0; row < n_rows; row++) {
			mpq_t coef;
			mpq_init(coef);
			mpq_QSget_coef(p_mpq, row, col, &coef);
			mpq_QSchange_coef(dual_p_mpq, col, row, coef);
			mpq_clear(coef);
		}
	}

	// Invert the objective sense for the dual problem
	int sense = 0;
	mpq_QSget_objsense(p_mpq, &sense);
	mpq_QSchange_objsense(dual_p_mpq, sense == QS_MIN ? QS_MAX : QS_MIN);

	for (int i = 0; i < n_cols; i++) {
		mpq_clear(obj_coeffs[i]);
	}
	for (int i = 0; i < n_rows; i++) {
		mpq_clear(rhss[i]);
	}
	mpq_clear(oneLpNum);

	return dual_p_mpq;
}

/* ========================================================================= */
/** @brief main function, here we build the LP, execute the options and exit */
int main (int argc,
					char **argv)
{
	int rval = 0;
	mpq_t v1,
	  v2;
	mpq_QSdata *p_mpq = 0;
	dbl_QSdata *p_dbl = 0;
	register unsigned i,
	  j,
	  k,
	  l;
	/* parse input */
	QSopt_ex_version();
	QSexactStart();
	rval = sl_parseargs (argc, argv);
	if (rval)
		return rval;
	mpq_init (v1);
	mpq_init (v2);
	/* create the problem with the appropriate number of variables and
	 * constraints */
	snprintf (strtmp, (size_t)1023, "OA_SL_%d-%d_%d-%d_t-%d", s1, k1, s2, k2, t);
	p_mpq = mpq_QScreate_prob (strtmp, QS_MIN);
	for (i = 0; i <= k1; i++)
		for (j = 0; j <= k2; j++)
		{
			snprintf (strtmp, (size_t)1023, "V_%d_%d", i, j);
			rval =
				mpq_QSnew_col (p_mpq, use_feas ? mpq_zeroLpNum : mpq_oneLpNum,
											 (i + j == 0) ? mpq_oneLpNum : mpq_zeroLpNum,
											 mpq_ILL_MAXDOUBLE, strtmp);
			CHECKRVALG (rval, CLEANUP);
			strtmp[0] = 'C';
			rval = mpq_QSnew_row (p_mpq, mpq_zeroLpNum,
											 ((i + j >= 1) && (i + j <= t)) ? 'E' : 'G', strtmp);
			CHECKRVALG (rval, CLEANUP);
		}
	/* now set the coefficients */
	for (i = 0; i <= k1; i++)
		for (j = 0; j <= k2; j++)
			for (k = 0; k <= k1; k++)
				for (l = 0; l <= k2; l++)
				{
					Kpoly (k, s1, i, k1, mpq_numref (v2));
					Kpoly (l, s2, j, k2, mpq_numref (v1));
					mpq_mul (v1, v1, v2);
					if (mpz_cmp_ui (mpq_numref (v1), (unsigned long)0))
					{
						rval =
							mpq_QSchange_coef (p_mpq, ((int)(j + i * (k2 + 1))), (int)(l + k * (k2 + 1)), v1);
						CHECKRVALG (rval, CLEANUP);
					}
				}
	/* now, if we are using scaling, we scale the non-zeros */
	if (use_scaling)
	{
		mpq_set_ui (v1, (unsigned long)1, (unsigned long)1);
		EXutilSimplify ((unsigned)(p_mpq->qslp->A.matsize), p_mpq->qslp->A.matval, v1);
		mpq_div (v2, mpq_oneLpNum, v1);
		fprintf (stderr, "Scale factor %lf\n", mpq_get_d (v2));
	}
	
	if (use_dual) {
		mpq_QSdata * dual_p_mpq = to_dual(p_mpq, strtmp);
		mpq_QSfree_prob(p_mpq);
		p_mpq = dual_p_mpq;
	}

	/* now we save the LP */
	if (use_double)
	{
		snprintf (strtmp, (size_t)1023, "OA_SL_%d-%d_%d-%d_t-%d", s1, k1, s2, k2, t);
		p_dbl = QScopy_prob_mpq_dbl (p_mpq, strtmp);
		rval = dbl_QSwrite_prob (p_dbl, out_file, use_mps ? "MPS" : "LP");
		CHECKRVALG (rval, CLEANUP);
	}
	else
	{
		rval = mpq_QSwrite_prob (p_mpq, out_file, use_mps ? "MPS" : "LP");
		CHECKRVALG (rval, CLEANUP);
	}

	/* ending */
CLEANUP:
	if (p_mpq)
		mpq_QSfree_prob (p_mpq);
	if (p_dbl)
		dbl_QSfree_prob (p_dbl);
	mpq_clear (v1);
	mpq_clear (v2);
	QSexactClear();
	return rval;
}
