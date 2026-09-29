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

#ifndef QS_DUMP_H__
#define QS_DUMP_H__

#include "qstruct_mpq.h"
#include "qstruct_rat.h"

void mpq_QSdump_xbz (const mpq_QSdata *p_mpq);
void mpq_QSdump_piz (const mpq_QSdata *p_mpq);
void mpq_QSdump_bz (const mpq_QSdata *p_mpq);
void mpq_QSdump_xnbz (const mpq_QSdata *p_mpq);
void mpq_QSdump_bfeas (const mpq_QSdata *p_mpq);
void mpq_QSdump_array (const mpq_t *array, const char* tag);

void mpq_QSdump_prob_col (const mpq_QSdata *p_mpq, int index, int col, char type);
void mpq_QSdump_prob (const mpq_QSdata *p_mpq);

// p_mpq can't be const because it is passed to mpq_ILLlib_tableau()
int mpq_QSdump_basis (mpq_QSdata *p_mpq);

void rat_QSdump_xbz (const rat_QSdata *p_rat);
void rat_QSdump_piz (const rat_QSdata *p_rat);
void rat_QSdump_bz (const rat_QSdata *p_rat);
void rat_QSdump_xnbz (const rat_QSdata *p_rat);
void rat_QSdump_bfeas (const rat_QSdata *p_rat);
void rat_QSdump_array (const rat_t *array, const char* tag);

void rat_QSdump_prob_col (const rat_QSdata *p_rat, int index, int col, char type);
void rat_QSdump_prob (const rat_QSdata *p_rat);

// p_rat can't be const because it is passed to rat_ILLlib_tableau()
int rat_QSdump_basis (rat_QSdata *p_rat);

#endif /* ! QS_DUMP_H__ */
