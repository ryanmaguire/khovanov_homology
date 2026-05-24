#ifndef HYBRIDKHOVANOV_H
#define HYBRIDKHOVANOV_H
#include "Komplex.h"
#include <stdbool.h>
/*
 * Compute the Khovanov chain complex using Algebraic reduction (Greedy Gaussian elimination).
 * 
 * If length <= threshold, uses Bar-Natan's iterative scan (tensoring one crossing and reducing).
 * If length > threshold, uses Divide-and-Conquer (recursively compute halves, tensor, reduce).
 */
Komplex *HybridKhovanov_compute_algebraic(int n_strands, int *crossings, bool *signs, int length, int threshold);
/*
 * Compute the Khovanov chain complex using Topological Morse Oracle reduction.
 * 
 * This builds the full unreduced chain complex for the entire braid, generates the
 * discrete Morse collapse schedule from the augmented rhomboid graph (Fixed-Strand),
 * and applies the schedule.
 */
Komplex *HybridKhovanov_compute_oracle(int n_strands, int *crossings, bool *signs, int length);
#endif // HYBRIDKHOVANOV_H
