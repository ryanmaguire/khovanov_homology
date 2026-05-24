#ifndef TANGLEKOMPLEX_H
#define TANGLEKOMPLEX_H
#include "Komplex.h"
#include "CobMatrix.h"
#include "CannedCobordismImpl.h"
/*
 * Create the unreduced Khovanov chain complex of a single crossing.
 * Positive crossing:
 *   0 -> [0-smoothing] -> [1-smoothing] -> 0
 * Negative crossing:
 *   0 -> [0-smoothing] -> [1-smoothing] -> 0
 * (with hpower and vpower shifts appropriate for Bar-Natan's normalization).
 */
Komplex *Komplex_singleCrossing(int n_strands, int crossing_index, bool positive);
/*
 * Create the length-1 Khovanov chain complex of the identity braid on n_strands.
 * It contains a single cap (the identity planar matching) in homological degree 0.
 */
Komplex *Komplex_identityBraid(int n_strands);
/*
 * Take the tensor product of two tangle chain complexes.
 * The two tangles are placed side-by-side and connected.
 * Their right endpoints (of k_left) are connected to the left endpoints (of k_right).
 */
Komplex *Komplex_compose_tangles(Komplex *k_left, Komplex *k_right, int n_strands);
#endif // TANGLEKOMPLEX_H
