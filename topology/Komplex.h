#ifndef KOMPLEX_H
#define KOMPLEX_H

#include <stdbool.h>
#include "CobMatrix.h"
#include "MorseOracle.h"

/*
 * Represents the Bar-Natan Khovanov Chain Complex
 *
 * C_0 --d_0--> C_1 --d_1--> C_2 --d_2--> ...
 *
 * differentials[i] maps C_i → C_{i+1}.
 * length = total number of chain groups.
 * differentials has size (length - 1).
 */
typedef struct Komplex {
    int length;               // Total number of chain columns
    CobMatrix **differentials; // Array of size length-1 mapping C_i -> C_{i+1}
} Komplex;

Komplex *Komplex_create(int length);
void Komplex_free(Komplex *k);

/*
 * Purpose:
 * Perform algebraic Gaussian elimination (Block Reduction Lemma)
 * on the differential at chain_idx.
 *
 * Given that D = differentials[chain_idx], and the entry
 * D[target_row][source_col] = φ is an isomorphism, this function:
 *
 * 1. Computes the Schur complement on D:
 *    D'[i][j] = D[i][j] - D[i][source_col] · φ^{-1} · D[target_row][j]
 *    for all i ≠ target_row, j ≠ source_col.
 *
 * 2. Removes target_row from D's target and source_col from D's source.
 *
 * 3. Updates the adjacent differential d_{chain_idx-1} (if exists):
 *    removes the row corresponding to source_col.
 *
 * 4. Updates the adjacent differential d_{chain_idx+1} (if exists):
 *    removes the column corresponding to target_row.
 *
 * Returns:
 * true if the reduction was performed successfully.
 */
bool Komplex_blockReductionLemma(Komplex *k, int chain_idx, int source_col, int target_row);

/*
 * Purpose:
 * Apply a full sequence of Gaussian eliminations guided by
 * the topological Morse oracle schedule.
 *
 * Arguments:
 * k — the Khovanov chain complex
 * schedule — the collapse schedule from MorseOracle_fullSchedule
 *
 * Returns:
 * Number of successful reductions applied.
 */
int Komplex_reduce_with_oracle(Komplex *k, const CollapseSchedule *schedule);

/*
 * Purpose:
 * Verify the chain complex condition: d² = 0.
 * For each consecutive pair of differentials, checks that
 * differentials[i+1] ∘ differentials[i] = 0.
 *
 * Returns:
 * true if d² = 0 holds everywhere.
 */
bool Komplex_verify_d_squared(const Komplex *k);

#endif // KOMPLEX_H
