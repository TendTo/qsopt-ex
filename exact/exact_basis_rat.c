/* ========================================================================= */
/* ESolver "Exact Mixed Integer Linear Solver" provides some basic structures
 * and algorithms commons in solving MIP's
 *
 * Copyright (C) 2005 Daniel Espinoza.
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
/* ========================================================================= */
/** @file
 * @ingroup Esolver */
/** @addtogroup Esolver */
/** @{ */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "exact_basis.h"
#include "dump.h"

#include <stdlib.h>
#include <string.h>

#include "logging-private.h"

#include "util.h"
#include "eg_timer.h"
#include "eg_exutil.h"
#include "except.h"

#include "basis_rat.h"
#include "editor_dbl.h"
#include "editor_mpf.h"
#include "fct_rat.h"
#include "lpdata_rat.h"
#include "simplex_rat.h"

static int load_and_factor_basis (rat_QSdata * p_rat,
                                  QSbasis * const basis,
                                  int *singular)
{
  int rval = 0;
  EGcallD(rat_QSload_basis (p_rat, basis));
  if (p_rat->cache)
  {
    rat_ILLlp_cache_free (p_rat->cache);
    rat_clear (p_rat->cache->val);
    ILL_IFFREE (p_rat->cache, rat_ILLlp_cache);
  }
  p_rat->qstatus = QS_LP_MODIFIED;
  if(p_rat->qslp->sinfo)
  {
    rat_ILLlp_sinfo_free(p_rat->qslp->sinfo);
    ILL_IFFREE(p_rat->qslp->sinfo, rat_ILLlp_sinfo);
  }
  if(p_rat->qslp->rA)
  {
    rat_ILLlp_rows_clear (p_rat->qslp->rA);
    ILL_IFFREE (p_rat->qslp->rA, rat_ILLlp_rows);
  }
  rat_free_internal_lpinfo (p_rat->lp);
  rat_init_internal_lpinfo (p_rat->lp);
  EGcallD(rat_build_internal_lpinfo (p_rat->lp));
  rat_ILLfct_set_variable_type (p_rat->lp);
  EGcallD(rat_ILLbasis_load (p_rat->lp, p_rat->basis));
  EGcallD(rat_ILLbasis_factor (p_rat->lp, singular));  // Always sets *singular
CLEANUP:
  return rval;
}

/* ========================================================================= */
/** @brief get the status for a given basis in rational arithmetic, it should
 * also leave everything set to get primal/dual solutions when needed.
 * */
int QSexact_basis_status_rat (rat_QSdata * p_rat,
                              int *status,
                              QSbasis * const basis,
                              const int msg_lvl,
                              int *const simplexalgo)
{
  int rval = 0,
      singular;
  rat_feas_info fi;
  EGtimer_t local_timer;
  rat_EGlpNumInitVar (fi.totinfeas);
  EGtimerReset (&local_timer);
  EGtimerStart (&local_timer);

  EGcallD(load_and_factor_basis (p_rat, basis, &singular));

  rat_ILLfct_compute_piz (p_rat->lp);
  rat_ILLfct_compute_dz (p_rat->lp);
  rat_ILLfct_compute_xbz (p_rat->lp);
  rat_ILLfct_check_pfeasible (p_rat->lp, &fi, rat_zeroLpNum);
  rat_ILLfct_check_dfeasible (p_rat->lp, &fi, rat_zeroLpNum);
  rat_ILLfct_set_status_values (p_rat->lp, fi.pstatus, fi.dstatus, PHASEII,
      PHASEII);
  if (p_rat->lp->basisstat.optimal)
  {
    *status = QS_LP_OPTIMAL;
    EGcallD(rat_QSgrab_cache (p_rat, QS_LP_OPTIMAL));
  }
  else if (p_rat->lp->basisstat.primal_infeasible
      || p_rat->lp->basisstat.dual_unbounded)
  {
    if (*status == QS_LP_INFEASIBLE && simplexalgo)
      *simplexalgo = PRIMAL_SIMPLEX;
    *status = QS_LP_INFEASIBLE;
    p_rat->lp->final_phase = PRIMAL_PHASEI;
    p_rat->lp->pIpiz = rat_EGlpNumAllocArray (p_rat->lp->nrows);
    rat_ILLfct_compute_phaseI_piz (p_rat->lp);
  }
  else if (p_rat->lp->basisstat.primal_unbounded)
    *status = QS_LP_UNBOUNDED;
  else
    *status = QS_LP_UNSOLVED;
  EGtimerStop (&local_timer);
  if(!msg_lvl)
  {
    MESSAGE(0, "Performing Rational Basic Solve on %s, %s, check"
        " done in %lg seconds, PS %s %lg, DS %s %lg", p_rat->name,
          *status == QS_LP_OPTIMAL    ? "RAT_optimal"
        : *status == QS_LP_INFEASIBLE ? "RAT_infeasible"
        : *status == QS_LP_UNBOUNDED  ? "RAT_unbounded"
                                      : "RAT_unsolved",
        local_timer.time,
          p_rat->lp->basisstat.primal_feasible   ? "F"
        : p_rat->lp->basisstat.primal_infeasible ? "I"
                                                 : "U",
          p_rat->lp->basisstat.primal_feasible   ? rat_to_d(p_rat->lp->objval)
        : p_rat->lp->basisstat.primal_infeasible ? rat_to_d(p_rat->lp->pinfeas)
                                                 : rat_to_d(p_rat->lp->objbound),
          p_rat->lp->basisstat.dual_feasible   ? "F"
        : p_rat->lp->basisstat.dual_infeasible ? "I"
                                               : "U",
          p_rat->lp->basisstat.dual_feasible   ? rat_to_d(p_rat->lp->dobjval)
        : p_rat->lp->basisstat.dual_infeasible ? rat_to_d(p_rat->lp->dinfeas)
                                               : rat_to_d(p_rat->lp->objbound));
  }
CLEANUP:
  rat_EGlpNumClearVar (fi.totinfeas);
  return rval;
}

/* ========================================================================= */
/** @brief get the status for a given basis in rational arithmetic, it should
 * also leave everything set to get primal/dual solutions when needed.
 * */
int QSdelta_basis_status_rat (rat_QSdata * p_rat,
                          int *status,
                          QSbasis * const basis,
                          const int msg_lvl,
                          int *const simplexalgo)
{
  int rval = 0,
      singular;
  rat_feas_info fi;
  EGtimer_t local_timer;
  rat_t zero;
  rat_EGlpNumInitVar (zero);  // Inits to zero
  rat_EGlpNumInitVar (fi.totinfeas);
  EGtimerReset (&local_timer);
  EGtimerStart (&local_timer);

  EGcallD(load_and_factor_basis (p_rat, basis, &singular));

  rat_ILLfct_compute_xbz (p_rat->lp);
  rat_ILLfct_check_pfeasible (p_rat->lp, &fi, rat_zeroLpNum);
  p_rat->lp->final_phase = PRIMAL_PHASEI;  // For rat_QSget_infeas_array
  p_rat->lp->pIpiz = rat_EGlpNumAllocArray (p_rat->lp->nrows);
  p_rat->lp->pIdz = rat_EGlpNumAllocArray (p_rat->lp->nnbasic);
  rat_ILLfct_compute_phaseI_piz (p_rat->lp);
  rat_ILLfct_compute_phaseI_dz (p_rat->lp);

  if (p_rat->simplex_display >= 2)
  {
    unsigned sz;
    if (p_rat->simplex_display >= 3)
    {
      rat_QSdump_prob(p_rat);
      rat_QSdump_basis(p_rat);
    }
    QSlog("QSdelta_basis_status: xnbz =");
    rat_QSdump_xnbz(p_rat);
    QSlog("QSdelta_basis_status: xbz =");
    rat_QSdump_xbz(p_rat);
    QSlog("QSdelta_basis_status: bfeas =");
    rat_QSdump_bfeas(p_rat);
    QSlog("QSdelta_basis_status: pIpiz =");
    rat_QSdump_array(p_rat->lp->pIpiz, "pIpiz");
    QSlog("QSdelta_basis_status: pIdz =");
    rat_QSdump_array(p_rat->lp->pIdz, "pIdz");
  }

  rat_ILLfct_check_pIdfeasible (p_rat->lp, &fi, zero);  // Requires zero to be non-const
  rat_ILLfct_set_status_values (p_rat->lp, fi.pstatus, fi.dstatus,
                                           PHASEII,    PHASEI);
  if (p_rat->lp->probstat.primal_feasible
   || p_rat->lp->probstat.primal_unbounded)
    *status = QS_LP_FEASIBLE;
  else if (p_rat->lp->probstat.primal_infeasible)
  {
    if (*status == QS_LP_INFEASIBLE && simplexalgo)
      *simplexalgo = PRIMAL_SIMPLEX;  // More efficient than dual, if infeas
    *status = QS_LP_INFEASIBLE;
  }
  else
    *status = QS_LP_UNSOLVED;
  EGtimerStop (&local_timer);
  if(!msg_lvl)
  {
    MESSAGE(0, "Performing Rational Basic Solve on %s, %s, check"
        " done in %lg seconds, PS %s %lg, DS %s %lg", p_rat->name,
          *status == QS_LP_FEASIBLE   ? "RAT_feasible"
        : *status == QS_LP_INFEASIBLE ? "RAT_infeasible"
                                      : "RAT_unsolved",
        local_timer.time,
          p_rat->lp->basisstat.primal_feasible   ? "F"
        : p_rat->lp->basisstat.primal_infeasible ? "I"
                                                 : "U",
          p_rat->lp->basisstat.primal_feasible   ? rat_to_d(p_rat->lp->objval)
        : p_rat->lp->basisstat.primal_infeasible ? rat_to_d(p_rat->lp->pinfeas)
                                                 : rat_to_d(p_rat->lp->objbound),
          p_rat->lp->basisstat.dual_feasible   ? "F"
        : p_rat->lp->basisstat.dual_infeasible ? "I"
                                               : "U",
          p_rat->lp->basisstat.dual_feasible   ? rat_to_d(p_rat->lp->dobjval)
        : p_rat->lp->basisstat.dual_infeasible ? rat_to_d(p_rat->lp->dinfeas)
                                               : rat_to_d(p_rat->lp->objbound));
  }
CLEANUP:
  rat_EGlpNumClearVar (fi.totinfeas);
  rat_EGlpNumClearVar (zero);
  return rval;
}

/* ========================================================================= */
/** @} */
/* end of exact_basis.c */

