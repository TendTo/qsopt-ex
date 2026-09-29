#ifndef __DELTA_RAT_H__
#define __DELTA_RAT_H__

#include "basicdefs.h"
#include "qstruct_rat.h"
#include "exact.h"

/** @file
 *  @ingroup DeltaSolver */
/** @defgroup DeltaSolver delta-complete solver
 *  @{ */

/* ========================================================================= */
/** @brief Given an rat_QSdata problem, solve the corresponding
 * delta-feasibility problem exactly.
 * @param p_mpq problem for which to determine delta-feasibility exactly.
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
int QSdelta_solver_rat (rat_QSdata * p_orig,
                    rat_t delta,
                    rat_t * const x,
                    rat_t * const y,
                    QSbasis * const ebasis,
                    unsigned * precision,
                    int simplexalgo,
                    int *status,
                    delta_callback_rat_t delta_callback,
                    void *callback_data);

/** @} */

/* ========================================================================= */
/** @brief copy the primal solution out of p_mpq.
 * @param x where to store the feasible primal solution.
 * @param p_mpq the problem data.
 * */
/* ========================================================================= */
int QSdelta_copy_x_rat (rat_t * const x, rat_t * y, const rat_QSdata * p_mpq);

#endif
