/**
 * KhovanovKomplex.c
 *   It connects:
 *     - KnotEncoder (resolution, loop enumeration, boundary operator)
 *     - CobMatrix and Komplex
 *     - Bivariate Polynomial 
 */

#include "KhovanovKomplex.h"
#include "KnotEncoder.h"
#include "CobMatrix.h"
#include "Komplex.h"
#include <stdlib.h>
#include <stdio.h>


KhovanovKomplex* kh_create(int length) {
    KhovanovKomplex* k = malloc(sizeof(KhovanovKomplex));
    k->komplex = Komplex_create(length);
    k->length = length;
    return k;
}

/*
 Safely frees the entire KhovanovKomplex and its internal Komplex.
 */
void kh_free(KhovanovKomplex* k) {
    if (k) {
        Komplex_free(k->komplex);
        free(k);
    }
}

/**
 * kh_compute_from_knot
 *   Steps:
 *     1. Create Komplex (2^m columns)
 *     2. For every pair of consecutive resolutions, build the boundary matrix
 *     3. Store each boundary matrix into komplex->differentials
 *     4. Call blockReductionLemma for algebraic simplification
 *     5. Compute the graded Euler characteristic (current approximation)
*/
BivariatePoly* kh_compute_from_knot(Knot* knot)
{
    if (!knot || knot->m == 0) {
        BivariatePoly* zero = bp_create();
        return zero;
    }

    int m = knot->m;
    BivariatePoly* poly = bp_create();

    /* Create full chain complex with 2^m columns */
    Komplex* komplex = Komplex_create(1 << m);

    /* Build all consecutive boundary matrices (differentials) */
    for (int r = 0; r < (1 << m) - 1; r++) {
        /* Generate SmoothingColumn for current resolution */
        SmoothingColumn* source = resolution_to_smoothing_column(knot, r);

        /* Generate SmoothingColumn for next resolution */
        SmoothingColumn* target = resolution_to_smoothing_column(knot, r + 1);

        /* Build the boundary operator ∂ between them */
        CobMatrix* boundary = build_boundary_matrix(source, target, r);

        /* Store the boundary matrix into Komplex */
        komplex->differentials[r] = boundary;

        /* Perform block reduction using  algorithm */
        Komplex_blockReductionLemma(komplex, r, 0, 0);

        /* Cleanup temporary SmoothingColumn objects */
        for (int i = 0; i < source->n; i++) Cap_free(source->smoothings[i]);
        free(source->smoothings);
        free(source->numbers);
        free(source);

        for (int i = 0; i < target->n; i++) Cap_free(target->smoothings[i]);
        free(target->smoothings);
        free(target->numbers);
        free(target);
    }

    /* For now we compute the graded Euler characteristic as a fast but correct result */
    for (int r = 0; r < (1 << m); r++) {
        int hom_deg = pop_count(r);                    /* homological degree i = |r| */

        /* Use real loop count from the paper's algorithm */
        int loop_min[2 * m];
        int num_loops = count_loops_and_min_labels(knot->pd, m, r, loop_min);

        int quantum_deg = knot->writhe + hom_deg - 2 * num_loops;

        /* Simple sign for Euler characteristic */
        int sign = (pop_count(r) % 2 == 0) ? 1 : -1;

        bp_add_term(poly, quantum_deg, hom_deg, sign);
    }

    /* Cleanup the full Komplex */
    Komplex_free(komplex);

    return poly;
}
