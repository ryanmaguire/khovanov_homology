#ifndef KHOVANOV_KOMPLEX_H
#define KHOVANOV_KOMPLEX_H

#include <stdbool.h>
#include "../encoder/KnotEncoder.h"
#include "../../topology/Komplex.h"

/**
 * @struct KhovanovKomplex
 * @brief Glue layer between topological knot resolutions and the algebraic
 *        chain complex (cube of resolutions → CobMatrix differentials).
 *
 * Note: the production Poincaré polynomial path is kh_poincare() in
 * polynomial/homology/, not this cube builder. This module remains available
 * for topology-side experiments and future deloop/Gauss integration.
 */
typedef struct {
    Knot *knot;                 /* topology source (no ownership; DO NOT FREE) */
    Komplex *komplex;

    /* Pre-computed / cached topology data */
    int *states_count;          /* states_count[h] = # resolutions at degree h */
    int **states_by_h;          /* states_by_h[h][i] = resolution bitmask */
    SmoothingColumn **columns;  /* columns[h] = Cap column for degree h */
} KhovanovKomplex;

KhovanovKomplex *KhovanovKomplex_alloc(Knot *knot);
void KhovanovKomplex_free(KhovanovKomplex *kh_complex);
bool KhovanovKomplex_build(KhovanovKomplex *kh_complex);

#endif /* KHOVANOV_KOMPLEX_H */
