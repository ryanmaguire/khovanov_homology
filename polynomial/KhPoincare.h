#ifndef KH_POINCARE_H
#define KH_POINCARE_H

#include <stdbool.h>
#include "../polynomial/BivariatePoly.h"
#include "../encoder/KnotEncoder.h"

/*
 * Enhanced-state integer Khovanov homology → Poincaré polynomial.
 *
 * Pipeline (productized from trefoil_integer_poincare.c):
 *   1. Enumerate resolutions r ∈ {0,1}^m and circle components
 *   2. Enhanced states = assignments of {1, X} to circles (bitmasks)
 *   3. Quantum grading (report convention, normalized):
 *        j = h + #circles − 2·#X + n_+ − 2·n_-
 *      Homological grading (normalized):
 *        r = h − n_-
 *      where h = popcount(r), n_+ = (m+w)/2, n_- = (m−w)/2
 *   4. Boundary maps from Frobenius merge/split + Jordan–Wigner sign
 *   5. Smith normal form per (h, j) block → image ranks + torsion
 *   6. β_{r,j} = dim C_{h,j} − rank(d_h) − rank(d_{h−1})
 *   7. Assemble Kh(q,t) = Σ t^r q^j β_{r,j} via bp_from_ranks
 *
 * Requirements on Knot:
 *   - pd[c][0..3] are 0-based edge labels in {0,...,2m−1}, each label twice
 *   - writhe is set correctly
 *   - m ≥ 0; practical limit m ≤ 12 (full cube 2^m)
 *
 * Does NOT implement scanning/delooping — that stays in topology/.
 */

/* Optional torsion report: list of (r, j, d) meaning a Z/d summand in H^{r,j}. */
typedef struct {
    int r; /* normalized homological degree */
    int j; /* quantum degree */
    int d; /* torsion order |d| > 1 */
} KhTorsionSummand;

typedef struct {
    KhTorsionSummand *items;
    int count;
    int capacity;
} KhTorsionReport;

void kh_torsion_report_init(KhTorsionReport *rep);
void kh_torsion_report_free(KhTorsionReport *rep);

/*
 * Compute the free Poincaré polynomial Kh(q,t).
 * Returns a newly allocated BivariatePoly (caller bp_free's).
 * On failure returns an empty polynomial (non-NULL if allocation works).
 */
BivariatePoly *kh_poincare(const Knot *knot);

/*
 * Same as kh_poincare, and if tors != NULL appends torsion summands
 * (r, j, d) with |d| > 1 discovered on Smith diagonals of the outgoing
 * differentials (placed on the target homological degree, matching the
 * original trefoil demo).
 */
BivariatePoly *kh_poincare_with_torsion(const Knot *knot, KhTorsionReport *tors);

/* Convenience: left-handed trefoil PD (0-based labels, writhe = -3). */
Knot *knot_create_left_trefoil(void);

#endif /* KH_POINCARE_H */
