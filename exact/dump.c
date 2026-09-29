/* dump - Functions for outputting internal data structures
 *
 * Copyright (C) 2020  Martin Sidaway
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include "dump.h"

#include "lib_mpq.h"
#include "lib_rat.h"
#include "logging-private.h"

#include <assert.h>

void mpq_QSdump_xbz (const mpq_QSdata *p_mpq)
{
  if (!p_mpq->lp->xbz)
  {
    QSlog ("xbz is unset");
    return;
  }
  assert (__EGlpNumArraySize (p_mpq->lp->xbz) == p_mpq->lp->nrows);
  for (int i = 0; i < p_mpq->lp->nrows; ++i)
  {
    QSlog ("%d: %g", p_mpq->lp->baz[i], mpq_get_d (p_mpq->lp->xbz[i]));
  }
}

void mpq_QSdump_piz (const mpq_QSdata *p_mpq)
{
  if (!p_mpq->lp->piz)
  {
    QSlog ("piz is unset");
    return;
  }
  assert (__EGlpNumArraySize (p_mpq->lp->piz) == p_mpq->lp->nrows);
  for (int i = 0; i < p_mpq->lp->nrows; ++i)
  {
    if (mpq_sgn (p_mpq->lp->piz[i]) != 0)
      QSlog ("%d: %g", i, mpq_get_d (p_mpq->lp->piz[i]));
  }
}

void mpq_QSdump_bz (const mpq_QSdata *p_mpq)
{
  if (!p_mpq->lp->bz)
  {
    QSlog ("bz is unset");
    return;
  }
  assert (__EGlpNumArraySize (p_mpq->lp->bz) >= p_mpq->lp->nrows);
  for (int i = 0; i < p_mpq->lp->nrows; ++i)
  {
    if (mpq_sgn (p_mpq->lp->bz[i]) != 0)
      QSlog ("%d: %g", i, mpq_get_d (p_mpq->lp->bz[i]));
  }
}

void mpq_QSdump_xnbz (const mpq_QSdata *p_mpq)
{
  //if (!p_mpq->lp->nbaz || !p_mpq->lp->vstat || !p_mpq->qslp->lower || !p_mpq->qslp->upper)
  if (!p_mpq->lp->nbaz || !p_mpq->lp->vstat || !p_mpq->lp->lz || !p_mpq->lp->uz)
  {
    QSlog ("Something needed to compute xnbz is unset");
    return;
  }
  for (int i = 0; i < p_mpq->lp->nnbasic; ++i)
  {
    int col = p_mpq->lp->nbaz[i];
    if (p_mpq->lp->vstat[col] == STAT_LOWER)
      QSlog ("%d: [%g", col, mpq_get_d (p_mpq->lp->lz[col]));
    else if (p_mpq->lp->vstat[col] == STAT_UPPER)
      QSlog ("%d: %g]", col, mpq_get_d (p_mpq->lp->uz[col]));
    else
      QSlog ("%d: 0", col);
  }
}

void mpq_QSdump_bfeas (const mpq_QSdata *p_mpq)
{
  if (!p_mpq->lp->bfeas)
  {
    QSlog ("bfeas is unset");
    return;
  }
  for (int i = 0; i < p_mpq->lp->nrows; ++i)
  {
    QSlog ("%d: %d", i, p_mpq->lp->bfeas[i]);
  }
}

void mpq_QSdump_array (const mpq_t *array, const char* tag)
{
  if (!array)
  {
    QSlog ("%s is unset", tag);
    return;
  }
  unsigned sz = __EGlpNumArraySize (array);
  for (int i = 0; i < sz; ++i)
  {
    QSlog ("%d: %g", i, mpq_get_d (array[i]));
  }
}

void mpq_QSdump_prob_col (const mpq_QSdata *p_mpq, int index, int col, char type)
{
    if (type == 'S')
    {
      //QSlog_nonl ("{%p, %p}: ", p_mpq->qslp->lower, p_mpq->qslp->upper);
      if (mpq_cmp (p_mpq->qslp->lower[col], mpq_NINFTY) <= 0)
        QSlog_nonl ("[-inf, ");
      else
        QSlog_nonl ("[%g, ", mpq_get_d (p_mpq->qslp->lower[col]));
      if (mpq_cmp (p_mpq->qslp->upper[col], mpq_INFTY) >= 0)
        QSlog_nonl ("inf]: ");
      else
        QSlog_nonl ("%g]: ", mpq_get_d (p_mpq->qslp->upper[col]));
    }
    else if (type == 'L')
    {
      // N.B. needs index rather than col
      QSlog_nonl ("(%c %g): ", p_mpq->qslp->sense[index], mpq_get_d (p_mpq->qslp->rhs[index]));
    }
    for (int j = p_mpq->qslp->A.matbeg[col];
         j < p_mpq->qslp->A.matbeg[col] + p_mpq->qslp->A.matcnt[col];
         ++j)
    {
      if (j > p_mpq->qslp->A.matbeg[col])
        QSlog_nonl (", ");
      QSlog_nonl ("%d=%g", p_mpq->qslp->A.matind[j],
                           mpq_get_d (p_mpq->qslp->A.matval[j]));
    }
    QSlog_nonl ("\n");
}

void mpq_QSdump_prob (const mpq_QSdata *p_mpq)
{
  QSlog ("mpq_QSdump_prob:");
  for (int i = 0; i < p_mpq->qslp->nstruct; ++i)
  {
    int col = p_mpq->qslp->structmap[i];
    QSlog_nonl ("struct col %d (%d): ", i, col);
    mpq_QSdump_prob_col (p_mpq, i, col, 'S');
  }
  for (int i = 0; i < p_mpq->qslp->nrows; ++i)
  {
    int col = p_mpq->qslp->rowmap[i];
    QSlog_nonl ("logical col %d (%d): ", i, col);
    mpq_QSdump_prob_col (p_mpq, i, col, 'L');
  }
}

int mpq_QSdump_basis (mpq_QSdata *p_mpq)
{
  QSlog ("mpq_QSdump_basis:");
  if (!p_mpq->lp || p_mpq->lp->basisid == -1)
  {
    QSlog ("No basis available");
    return E_GENERAL_ERROR;
  }
  mpq_t *row = mpq_EGlpNumAllocArray (p_mpq->qslp->nstruct + p_mpq->qslp->nrows);
  int rval = 0;
  for (int i = 0; i < p_mpq->qslp->nrows; ++i)
  {
    EGcallD (mpq_ILLlib_tableau (p_mpq->lp, i, 0, row));
    QSlog_nonl ("row %d:", i);
    for (int j = 0; j < p_mpq->qslp->nstruct; ++j)
    {
      QSlog_nonl (" %g", mpq_get_d (row[j]));
    }
    QSlog_nonl (" |");
    for (int j = 0; j < p_mpq->qslp->nrows; ++j)
    {
      QSlog_nonl (" %g", mpq_get_d (row[p_mpq->qslp->nstruct + j]));
    }
    QSlog_nonl ("\n");
  }
CLEANUP:
  mpq_EGlpNumFreeArray (row);
  return rval;
}

void rat_QSdump_xbz (const rat_QSdata *p_rat)
{
  if (!p_rat->lp->xbz)
  {
    QSlog ("xbz is unset");
    return;
  }
  assert (__EGlpNumArraySize (p_rat->lp->xbz) == p_rat->lp->nrows);
  for (int i = 0; i < p_rat->lp->nrows; ++i)
  {
    QSlog ("%d: %g", p_rat->lp->baz[i], rat_to_d (p_rat->lp->xbz[i]));
  }
}

void rat_QSdump_piz (const rat_QSdata *p_rat)
{
  if (!p_rat->lp->piz)
  {
    QSlog ("piz is unset");
    return;
  }
  assert (__EGlpNumArraySize (p_rat->lp->piz) == p_rat->lp->nrows);
  for (int i = 0; i < p_rat->lp->nrows; ++i)
  {
    if (rat_sgn (p_rat->lp->piz[i]) != 0)
      QSlog ("%d: %g", i, rat_to_d (p_rat->lp->piz[i]));
  }
}

void rat_QSdump_bz (const rat_QSdata *p_rat)
{
  if (!p_rat->lp->bz)
  {
    QSlog ("bz is unset");
    return;
  }
  assert (__EGlpNumArraySize (p_rat->lp->bz) >= p_rat->lp->nrows);
  for (int i = 0; i < p_rat->lp->nrows; ++i)
  {
    if (rat_sgn (p_rat->lp->bz[i]) != 0)
      QSlog ("%d: %g", i, rat_to_d (p_rat->lp->bz[i]));
  }
}

void rat_QSdump_xnbz (const rat_QSdata *p_rat)
{
  //if (!p_rat->lp->nbaz || !p_rat->lp->vstat || !p_rat->qslp->lower || !p_rat->qslp->upper)
  if (!p_rat->lp->nbaz || !p_rat->lp->vstat || !p_rat->lp->lz || !p_rat->lp->uz)
  {
    QSlog ("Something needed to compute xnbz is unset");
    return;
  }
  for (int i = 0; i < p_rat->lp->nnbasic; ++i)
  {
    int col = p_rat->lp->nbaz[i];
    if (p_rat->lp->vstat[col] == STAT_LOWER)
      QSlog ("%d: [%g", col, rat_to_d (p_rat->lp->lz[col]));
    else if (p_rat->lp->vstat[col] == STAT_UPPER)
      QSlog ("%d: %g]", col, rat_to_d (p_rat->lp->uz[col]));
    else
      QSlog ("%d: 0", col);
  }
}

void rat_QSdump_bfeas (const rat_QSdata *p_rat)
{
  if (!p_rat->lp->bfeas)
  {
    QSlog ("bfeas is unset");
    return;
  }
  for (int i = 0; i < p_rat->lp->nrows; ++i)
  {
    QSlog ("%d: %d", i, p_rat->lp->bfeas[i]);
  }
}

void rat_QSdump_array (const rat_t *array, const char* tag)
{
  if (!array)
  {
    QSlog ("%s is unset", tag);
    return;
  }
  unsigned sz = __EGlpNumArraySize (array);
  for (int i = 0; i < sz; ++i)
  {
    QSlog ("%d: %g", i, rat_to_d (array[i]));
  }
}

void rat_QSdump_prob_col (const rat_QSdata *p_rat, int index, int col, char type)
{
    if (type == 'S')
    {
      //QSlog_nonl ("{%p, %p}: ", p_rat->qslp->lower, p_rat->qslp->upper);
      if (rat_cmp (p_rat->qslp->lower[col], rat_NINFTY) <= 0)
        QSlog_nonl ("[-inf, ");
      else
        QSlog_nonl ("[%g, ", rat_to_d (p_rat->qslp->lower[col]));
      if (rat_cmp (p_rat->qslp->upper[col], rat_INFTY) >= 0)
        QSlog_nonl ("inf]: ");
      else
        QSlog_nonl ("%g]: ", rat_to_d (p_rat->qslp->upper[col]));
    }
    else if (type == 'L')
    {
      // N.B. needs index rather than col
      QSlog_nonl ("(%c %g): ", p_rat->qslp->sense[index], rat_to_d (p_rat->qslp->rhs[index]));
    }
    for (int j = p_rat->qslp->A.matbeg[col];
         j < p_rat->qslp->A.matbeg[col] + p_rat->qslp->A.matcnt[col];
         ++j)
    {
      if (j > p_rat->qslp->A.matbeg[col])
        QSlog_nonl (", ");
      QSlog_nonl ("%d=%g", p_rat->qslp->A.matind[j],
                           rat_to_d (p_rat->qslp->A.matval[j]));
    }
    QSlog_nonl ("\n");
}

void rat_QSdump_prob (const rat_QSdata *p_rat)
{
  QSlog ("rat_QSdump_prob:");
  for (int i = 0; i < p_rat->qslp->nstruct; ++i)
  {
    int col = p_rat->qslp->structmap[i];
    QSlog_nonl ("struct col %d (%d): ", i, col);
    rat_QSdump_prob_col (p_rat, i, col, 'S');
  }
  for (int i = 0; i < p_rat->qslp->nrows; ++i)
  {
    int col = p_rat->qslp->rowmap[i];
    QSlog_nonl ("logical col %d (%d): ", i, col);
    rat_QSdump_prob_col (p_rat, i, col, 'L');
  }
}

int rat_QSdump_basis (rat_QSdata *p_rat)
{
  QSlog ("rat_QSdump_basis:");
  if (!p_rat->lp || p_rat->lp->basisid == -1)
  {
    QSlog ("No basis available");
    return E_GENERAL_ERROR;
  }
  rat_t *row = rat_EGlpNumAllocArray (p_rat->qslp->nstruct + p_rat->qslp->nrows);
  int rval = 0;
  for (int i = 0; i < p_rat->qslp->nrows; ++i)
  {
    EGcallD (rat_ILLlib_tableau (p_rat->lp, i, 0, row));
    QSlog_nonl ("row %d:", i);
    for (int j = 0; j < p_rat->qslp->nstruct; ++j)
    {
      QSlog_nonl (" %g", rat_to_d (row[j]));
    }
    QSlog_nonl (" |");
    for (int j = 0; j < p_rat->qslp->nrows; ++j)
    {
      QSlog_nonl (" %g", rat_to_d (row[p_rat->qslp->nstruct + j]));
    }
    QSlog_nonl ("\n");
  }
CLEANUP:
  rat_EGlpNumFreeArray (row);
  return rval;
}
