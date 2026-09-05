/*
 * Thin regression driver for step 3: left trefoil via kh_poincare().
 * Compile (from repo root, after placing files):
 *   gcc -O2 -o test_kh_poincare \
 *     polynomial/homology/KhPoincare.c \
 *     polynomial/polynomial/BivariatePoly.c \
 *     IntegerMatrix.c \
 *     polynomial/test/test_kh_poincare_trefoil.c \
 *     -I. -Ipolynomial/homology -Ipolynomial/polynomial -Ipolynomial/encoder
 *
 * Note: KnotEncoder.c is NOT required; Knot struct is in the header only.
 */

#include "KhPoincare.h"
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    Knot *knot = knot_create_left_trefoil();
    if (!knot) {
        fprintf(stderr, "failed to create left trefoil\n");
        return 1;
    }

    KhTorsionReport tors;
    kh_torsion_report_init(&tors);

    BivariatePoly *poly = kh_poincare_with_torsion(knot, &tors);
    printf("left trefoil Kh(q,t) = ");
    bp_print(poly);

    if (tors.count > 0) {
        printf("torsion summands:\n");
        for (int i = 0; i < tors.count; i++) {
            printf("  Z_%d at (r=%d, j=%d)\n",
                   tors.items[i].d, tors.items[i].r, tors.items[i].j);
        }
    } else {
        printf("torsion summands: (none recorded)\n");
    }

    BivariatePoly *chi = bp_eval_t(poly, -1);
    printf("chi_q (t=-1) = ");
    bp_print(chi);

    bp_free(chi);
    bp_free(poly);
    kh_torsion_report_free(&tors);
    knot_free(knot);
    return 0;
}
