#ifndef KOMPLEX_H
#define KOMPLEX_H

#include "CobMatrix.h"

/*
 * Represents the Bar-Natan Khovanov Chain Complex
 */
typedef struct Komplex {
  int length;                // Total number of chain columns
  CobMatrix **differentials; // Array of size length-1 mapping C_i -> C_{i+1}
} Komplex;

Komplex *Komplex_create(int length);
void Komplex_free(Komplex *k);

/*
 * Perform algebraic Gaussian elimination locally between chain i and i+1
 * guided by the deterministic topological Morse collapse schedule.
 * Returns true if successful reduction occurs.
 */
bool Komplex_blockReductionLemma(Komplex *k, int chain_idx, int source_col,
                                 int target_row);

#endif // KOMPLEX_H
