/** @file
 *  @ingroup DeltaSolver */
/** @addtogroup DeltaSolver
 *  @{ */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "delta_rat.h"
#include "exact_basis_rat.h"

#include "qstruct_dbl.h"
#include "qstruct_mpf.h"

#include "basis_rat.h"
#include "except.h"
#include "editor_dbl.h"
#include "editor_mpf.h"
#include "eg_macros.h"
#include "eg_timer.h"
#include "fct_rat.h"
#include "qsopt_mpf.h"
#include "qsopt_dbl.h"
#include "simplex_rat.h"
#include "dump.h"


int QSdelta_copy_x_rat (rat_t * const x, rat_t * y, const rat_QSdata * p_rat)
{
  int rval = 0;
  int i, col;
  rat_t *tempx = 0;
  rat_lpinfo* lp = p_rat->lp;
  rat_ILLlpdata* qslp = p_rat->qslp;

  // Populate x with values of structural variables
  if (lp->nrows != qslp->nrows ||
      lp->ncols != qslp->ncols ||
      lp->nnbasic != qslp->nstruct ||
      lp->ncols != lp->nrows + lp->nnbasic)
  {
    QSlog("Unexpected condition: lp and qslp dimensions do not match");
    rval = 1;
    ILL_CLEANUP;
  }
  tempx = rat_EGlpNumAllocArray (lp->ncols);
  // Set basic variables
  for (i = 0; i < lp->nrows; i++)
    rat_EGlpNumCopy (tempx[lp->baz[i]], lp->xbz[i]);
  // Set non-basic variables
  for (i = 0; i < lp->nnbasic; i++)
  {
    col = lp->nbaz[i];
    if (lp->vstat[col] == STAT_UPPER)
      rat_EGlpNumCopy (tempx[col], lp->uz[col]);
    else if (lp->vstat[col] == STAT_LOWER)
      rat_EGlpNumCopy (tempx[col], lp->lz[col]);
    else
      rat_EGlpNumZero (tempx[col]);
  }
  // Get structural variables
  for (i = 0; i < qslp->nstruct; i++)
  {
    rat_EGlpNumCopy (x[i], tempx[qslp->structmap[i]]);
  }
  if (y && lp && lp->pIpiz)
  {
    for (i = 0; i < qslp->nrows; i++)
    {
      rat_EGlpNumCopy (y[i], lp->pIpiz[i]);
    }
  }

CLEANUP:

  rat_EGlpNumFreeArray (tempx);

  EG_RETURN (rval);
}

/* ========================================================================= */
/** @brief print into screen (if enable) a message indicating that we have
 * successfully proven infeasibility.
 * @param p_rat the problem data.
 * */
/* ========================================================================= */
static int infeasible_output (rat_QSdata * p_rat,
                              rat_t * const x, 
                             rat_t * const y)
{
  int rval = 0;

  if (p_rat->simplex_display)
    QSlog("Problem is infeasible");
  if (x)
    EGcallD(QSdelta_copy_x_rat (x, y, p_rat));

CLEANUP:

  EG_RETURN (rval);
}

/* ========================================================================= */
/** @brief print into screen (if enable) a message indicating that we have
 * successfully proven feasibility.
 * @param p_rat the problem data.
 * */
/* ========================================================================= */
static int feasible_output (rat_QSdata * p_rat,
                            rat_t * const x, 
                            rat_t * const y)
{
  int rval = 0;

  if (p_rat->simplex_display)
    QSlog("Problem is feasible");
  if (x)
    EGcallD(QSdelta_copy_x_rat (x, y, p_rat));

CLEANUP:

  EG_RETURN (rval);
}

/* ========================================================================= */
/** Check for and report delta-feasibility of the basis
 * @param p_rat the problem data.
 * @param delta the maximum infeasibility of a delta-feasible solution; updated
 * with the actual infeasibility if delta-feasibility determined.
 * @param status where to store the new status (feasible or delta-feasible),
 * if applicable
 * @param x where to store a delta-feasible primal solution (if not null).
 * */
/* ========================================================================= */
static int check_delta_feas (rat_QSdata const * p_rat,
                             rat_t delta,
                             int *status,
                             rat_t * const x,
                             rat_t * const y,
                             delta_callback_rat_t delta_callback,
                             rat_t last_infeas,
                             void *callback_data)
{
  int i, col;
  rat_t infeas, err1, err2;
  int rval = 0;
  rat_lpinfo* lp = p_rat->lp;
  rat_ILLlpdata* qslp = p_rat->qslp;

  *status = QS_LP_UNSOLVED;

  rat_EGlpNumInitVar (infeas);
  rat_EGlpNumInitVar (err1);
  rat_EGlpNumInitVar (err2);
  rat_EGlpNumZero (infeas);

  for (i = 0; i < lp->nrows; i++)
  {
    col = lp->baz[i];
    rat_EGlpNumCopyDiff (err1, lp->xbz[i], lp->uz[col]);
    rat_EGlpNumCopyDiff (err2, lp->lz[col], lp->xbz[i]);
    if (rat_EGlpNumIsLess (rat_zeroLpNum, err1)
        && rat_EGlpNumIsNeq (lp->uz[col], rat_INFTY, rat_zeroLpNum))
    {
      if (rat_EGlpNumIsLess (infeas, err1))
        rat_EGlpNumCopy (infeas, err1);
      WARNINGL (QSE_WLVL, rat_EGlpNumIsLess (rat_INFTY, err1),
               "This is impossible: lu = %15lg xbz = %15lg" " rat_INFTY = %15lg",
               rat_EGlpNumToLf (lp->uz[col]), rat_EGlpNumToLf (lp->xbz[i]),
               rat_EGlpNumToLf (rat_INFTY));
    }
    else if (rat_EGlpNumIsLess (rat_zeroLpNum, err2)
             && rat_EGlpNumIsNeq (lp->lz[col], rat_NINFTY, rat_zeroLpNum))
    {
      if (rat_EGlpNumIsLess (infeas, err2))
        rat_EGlpNumAddTo (infeas, err2);
      WARNINGL (QSE_WLVL, rat_EGlpNumIsLess (rat_INFTY, err2),
               "This is impossible: lz = %15lg xbz = %15lg" " rat_NINFTY = %15lg",
               rat_EGlpNumToLf (lp->lz[col]), rat_EGlpNumToLf (lp->xbz[i]),
               rat_EGlpNumToLf (rat_NINFTY));
    }
  }

  if (rat_EGlpNumIsLessZero (infeas))
  {
    QSlog("Negative infeasibility (impossible): %lf %la",
                rat_EGlpNumToLf (infeas), rat_EGlpNumToLf (infeas));
  }

  if (!rat_EGlpNumIsNeqqZero (infeas))
  {
    // feasible
    if (p_rat->simplex_display)
    {
      QSlog("Problem is feasible");
    }
    rat_EGlpNumCopy (delta, infeas);
    *status = QS_LP_FEASIBLE;
  }
  else if (!rat_EGlpNumIsLess (delta, infeas))
  {
    // delta-feasible
    if (p_rat->simplex_display)
    {
      QSlog("Problem is delta-feasible with delta = %lf",
            rat_EGlpNumToLf (infeas));
    }
    rat_EGlpNumCopy (delta, infeas);
    *status = QS_LP_DELTA_FEASIBLE;
  }
  else if (NULL != delta_callback)
  {
    if (rat_sgn (last_infeas) == 0 || rat_cmp (infeas, last_infeas) < 0)
    {
      rat_set (last_infeas, infeas);
      delta_callback(p_rat, x, infeas, delta, callback_data);
    }
  }

  if (x && (QS_LP_FEASIBLE == *status || QS_LP_DELTA_FEASIBLE == *status))
    EGcallD(QSdelta_copy_x_rat (x, y, p_rat));

  if (QS_LP_FEASIBLE != *status && QS_LP_DELTA_FEASIBLE != *status)
  {
    IFMESSAGE(p_rat->simplex_display, "Failed to make final conclusion on basis of solver result");
  }

CLEANUP:

  rat_EGlpNumClearVar (infeas);
  rat_EGlpNumClearVar (err1);
  rat_EGlpNumClearVar (err2);

  EG_RETURN (rval);
}

/* ========================================================================= */
/** @brief Given an rat_QSdata problem, solve the corresponding
 * delta-feasibility problem exactly.
 * @param p_rat problem for which to determine delta-feasibility exactly.
 * @param delta the delta to use for determining delta-feasibility; the maximum
 * perturbation of RHS/bounds required to make a delta-feasible solution
 * feasible; updated with the actual infeasibility if delta-feasibility
 * determined.
 * @param x if not null, we store here a delta-feasible solution to the problem
 * (if delta-feasibility established).
 * @param y if not null, we store here a certificate of infeasibility for the
 * problem (if infeasibility established).
 * @param ebasis if not null, use the given basis to start the iteration of
 * simplex, and store here the final basis (where applicable).
 * @param precision pointer to the variable where we will store the actual precision used to solve the problem; 
 * can be null, in which case we will not report the actual precision used.
 * @param simplexalgo whether to use primal or dual simplex while solving the
 * delta-feasibility problem.
 * @param status pointer to the integer where we will return the status of the
 * problem, either feasible, delta-feasible or infeasible (we could also return
 * time out).
 * @param delta_callback if not null, will be called if a delta-satisfying
 * result is found for some value greater than delta.
 * @param callback_data additional parameter to be passed to delta_callback.
 * @return zero on success, non-zero otherwise. */
int QSdelta_solver_rat (rat_QSdata * p_rat,
                    rat_t delta,
                    rat_t * const x,
                    rat_t * const y,
                    QSbasis * const ebasis,
                    unsigned * precision,
                    int simplexalgo,
                    int *status,
                    delta_callback_rat_t delta_callback,
                    void *callback_data)
{
  /* local variables */
  int last_status = 0, last_iter = 0;
  QSbasis *basis = 0;
  unsigned fallback_precision;
  if (!precision) precision = &fallback_precision;
  *precision = sizeof(double) * 8;
  int rval = 0,
    it = QS_EXACT_MAX_ITER;
  dbl_QSdata *p_dbl = 0;
  mpf_QSdata *p_mpf = 0;
  double *x_dbl = 0,
   *y_dbl = 0;
  mpf_t *x_mpf = 0,
   *y_mpf = 0;
  rat_t last_infeas;
  rat_init (last_infeas);
  int const msg_lvl = __QS_SB_VERB <= DEBUG ? 0: (1 - p_rat->simplex_display) * 10000;
  *status = QS_LP_UNSOLVED;
  /* set the objective function to zero (in the copy) (enforced in the preconditions) */
  // EGcallD (rat_QSclear_obj (p_rat));
  /* save the problem if we are really debugging */
  if(DEBUG >= __QS_SB_VERB)
  {
    EGcallD(rat_QSwrite_prob(p_rat, "qsxprob.lp","LP"));
  }
  /* try first with doubles */
  if (p_rat->simplex_display || DEBUG >= __QS_SB_VERB)
  {
    QSlog("Trying double precision");
  }
  p_dbl = QScopy_prob_rat_dbl (p_rat, "dbl_problem");
  if(__QS_SB_VERB <= DEBUG && !p_dbl->simplex_display) p_dbl->simplex_display = 1;
  if (ebasis && ebasis->nstruct)
    dbl_QSload_basis (p_dbl, ebasis);
  if (dbl_ILLeditor_solve (p_dbl, simplexalgo))
  {
    MESSAGE(p_rat->simplex_display ? 0: __QS_SB_VERB,
            "double approximation failed, code %d, "
            "continuing in extended precision", rval);
    goto MPF_PRECISION;
  }
  EGcallD(dbl_QSget_status (p_dbl, status));
  if ((*status == QS_LP_INFEASIBLE) &&
      (p_dbl->lp->final_phase != PRIMAL_PHASEI) &&
      (p_dbl->lp->final_phase != DUAL_PHASEII))
    dbl_QSopt_primal (p_dbl, status);
  EGcallD(dbl_QSget_status (p_dbl, status));
  last_status = *status;
  EGcallD(dbl_QSget_itcnt(p_dbl, 0, 0, 0, 0, &last_iter));
  /* deal with the problem depending on what status we got from our optimizer */
  if (QS_LP_OPTIMAL == *status || QS_LP_UNBOUNDED == *status || QS_LP_INFEASIBLE == *status)
  {
    basis = dbl_QSget_basis (p_dbl);
    MESSAGE (msg_lvl, "Basis hash is 0x%016lX", QSexact_basis_hash(basis));
    EGcallD(QSdelta_basis_status_rat (p_rat, status, basis, msg_lvl, &simplexalgo));
    if (QS_LP_INFEASIBLE == *status)
    {
      infeasible_output (p_rat, x, y);
      goto CLEANUP;
    }
    else if (QS_LP_FEASIBLE == *status)
    {
      EGcallD(feasible_output (p_rat, x, y));
      rat_EGlpNumCopy (delta, rat_zeroLpNum);
      goto CLEANUP;
    }
    /* check for delta-feasibility */
    EGcallD(check_delta_feas (p_rat, delta, status, x, y, delta_callback, last_infeas, callback_data));
    if (QS_LP_FEASIBLE == *status || QS_LP_DELTA_FEASIBLE == *status)
      goto CLEANUP;
  }
  else
  {
    MESSAGE (msg_lvl, "Floating-point solver reports failure; trying next precision");
  }
  IFMESSAGE(p_rat->simplex_display,"Retrying in extended precision");
  /* if we reach this point, then we have to keep going, we use the previous
   * basis ONLY if the previous precision thinks that it has the optimal
   * solution, otherwise we start from scratch. */
  *precision = 128;
  MPF_PRECISION:
  dbl_QSfree_prob (p_dbl);
  p_dbl = 0;
  /* try with multiple precision floating points */
  for (; it--; *precision = (unsigned) (*precision * 1.5))
  {
    QSexact_set_precision (*precision);
    if (p_rat->simplex_display || DEBUG >= __QS_SB_VERB)
    {
      QSlog("Trying mpf with %u bits", *precision);
    }
    p_mpf = QScopy_prob_rat_mpf (p_rat, "mpf_problem");
    if(DEBUG >= __QS_SB_VERB)
    {
      EGcallD(mpf_QSwrite_prob(p_mpf, "qsxprob.mpf.lp","LP"));
    }
    if(__QS_SB_VERB <= DEBUG && !p_mpf->simplex_display) p_mpf->simplex_display = 1;
    simplexalgo = PRIMAL_SIMPLEX;
    if(!last_iter) last_status = QS_LP_UNSOLVED;
    if(last_status == QS_LP_OPTIMAL || last_status == QS_LP_INFEASIBLE)
    {
      if (p_rat->simplex_display || DEBUG >= __QS_SB_VERB)
      {
        QSlog("Reusing previous basis");
      }
      if (basis)
      {
        EGcallD(mpf_QSload_basis (p_mpf, basis));
        mpf_QSfree_basis (basis);
        simplexalgo = PRIMAL_SIMPLEX;
        basis = 0;
      }
      else if (ebasis && ebasis->nstruct)
      {
        mpf_QSload_basis (p_mpf, ebasis);
        simplexalgo = PRIMAL_SIMPLEX;
      }
    }
    else
    {
      if(p_mpf->basis)
      {
        mpf_ILLlp_basis_free(p_mpf->basis);
        p_mpf->lp->basisid = -1;
        p_mpf->factorok = 0;
      }
      if (p_rat->simplex_display || DEBUG >= __QS_SB_VERB)
      {
        QSlog("Not using previous basis");
      }
    }
    if (mpf_ILLeditor_solve (p_mpf, simplexalgo))
    {
      if (p_rat->simplex_display || DEBUG >= __QS_SB_VERB)
      {
        QSlog("mpf_%u precision failed, error code %d, continuing with "
                    "next precision", *precision, rval);
       }
      goto NEXT_PRECISION;
    }
    EGcallD(mpf_QSget_status (p_mpf, status));
    if ((*status == QS_LP_INFEASIBLE) &&
        (p_mpf->lp->final_phase != PRIMAL_PHASEI) &&
        (p_mpf->lp->final_phase != DUAL_PHASEII))
      mpf_QSopt_primal (p_mpf, status);
    EGcallD(mpf_QSget_status (p_mpf, status));
    last_status = *status;
    EGcallD(mpf_QSget_itcnt(p_mpf, 0, 0, 0, 0, &last_iter));
    /* deal with the problem depending on status we got from our optimizer */
    if (QS_LP_OPTIMAL == *status || QS_LP_UNBOUNDED == *status || QS_LP_INFEASIBLE == *status)
    {
      mpf_QSfree_basis (basis);
      basis = mpf_QSget_basis (p_mpf);
      MESSAGE (msg_lvl, "Basis hash is 0x%016lX", QSexact_basis_hash(basis));
      EGcallD(QSdelta_basis_status_rat (p_rat, status, basis, msg_lvl, &simplexalgo));
      if (QS_LP_INFEASIBLE == *status)
      {
        infeasible_output (p_rat, x, y);
        goto CLEANUP;
      }
      else if (QS_LP_FEASIBLE == *status)
      {
        EGcallD(feasible_output (p_rat, x, y));
        rat_EGlpNumCopy (delta, rat_zeroLpNum);
        goto CLEANUP;
      }
      /* check for delta-feasibility */
      EGcallD(check_delta_feas (p_rat, delta, status, x, y, delta_callback, last_infeas, callback_data));
      if (QS_LP_FEASIBLE == *status || QS_LP_DELTA_FEASIBLE == *status)
        goto CLEANUP;
    }
    else
    {
      MESSAGE (msg_lvl, "Floating-point solver reports failure; trying next precision");
    }
  NEXT_PRECISION:
    mpf_QSfree_prob (p_mpf);
    p_mpf = 0;
  }
  /* ending */
CLEANUP:
  rat_clear (last_infeas);
  dbl_EGlpNumFreeArray (x_dbl);
  dbl_EGlpNumFreeArray (y_dbl);
  mpf_EGlpNumFreeArray (x_mpf);
  mpf_EGlpNumFreeArray (y_mpf);
  if (ebasis && basis)
  {
    ILL_IFFREE (ebasis->cstat, char);
    ILL_IFFREE (ebasis->rstat, char);
    ebasis->nstruct = basis->nstruct;
    ebasis->nrows = basis->nrows;
    ebasis->cstat = basis->cstat;
    ebasis->rstat = basis->rstat;
    basis->cstat = basis->rstat = 0;
  }
  rat_QSfree_basis (basis);
  dbl_QSfree_prob (p_dbl);
  mpf_QSfree_prob (p_mpf);

  EG_RETURN (rval);
}


/** @} */
