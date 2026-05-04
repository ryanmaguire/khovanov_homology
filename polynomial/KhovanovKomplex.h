#ifndef KHOVANOV_KOMPLEX_H
#define KHOVANOV_KOMPLEX_H

#include "IntMatrix.h"
#include "KnotEncoder.h"
#include "BivariatePoly.h"

// Full Khovanov chain complex using IntMatrix differentials
typedef struct {
    int ncolumns;           // homological degrees 0 to m
    IntMatrix** differentials;  // d_i : C_i → C_{i-1}
    int startnum;           // overall shift
} KhovanovKomplex;

// Create complex from knot PD and build full chain + Smith form
KhovanovKomplex* kh_create_from_knot(Knot* knot);

// Compute Poincaré polynomial (graded Betti numbers)
BivariatePoly* kh_compute_poincare(KhovanovKomplex* k);

// Free everything
void kh_free(KhovanovKomplex* k);

#endif