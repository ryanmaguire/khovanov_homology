#ifndef TOPOLOGY_BRIDGE_H
#define TOPOLOGY_BRIDGE_H

/*
 * Topology ↔ polynomial interface contract
 * ========================================
 *
 * OWNERSHIP BOUNDARY
 * ------------------
 * topology/     owns: tangle scanning, delooping, Gaussian elimination on
 *               CobMatrix / CannedCobordism, evaluation of closed cobordisms
 *               to integer coefficients, and production of bigraded integer
 *               boundary matrices (or already-reduced rank tables).
 *
 * polynomial/   owns: Smith extraction (via IntegerMatrix), Betti assembly,
 *               BivariatePoly arithmetic, Poincaré output, Jones self-check.
 *               MUST NOT reimplement scanning / deloop / cobordism composition.
 *
 * CURRENT PRODUCTION PATH (steps 1–5)
 * -----------------------------------
 *   Knot (0-based PD)
 *     → kh_poincare()          [polynomial/homology/KhPoincare.c]
 *     → enhanced states + full cube + IntegerMatrix Smith
 *     → BivariatePoly Kh(q,t)
 *
 * This is correct for small m (≤ ~12–16) but is NOT FastKh.
 *
 * TARGET FAST PATH (topology delivers, polynomial consumes)
 * ---------------------------------------------------------
 *   Knot / braid word
 *     → topology FullScanning / HybridKhovanov
 *         (compose local complexes → deloop → Gauss elim)
 *     → for each quantum degree j (and each homological step h):
 *           integer Mat *d_h^{(j)}   already in Z coefficients
 *        OR precomputed image ranks + torsion of those matrices
 *     → polynomial:
 *           smith_extract / toSmithForm_extract
 *           β_{r,j} = dim C_{h,j} − rank(d_h) − rank(d_{h−1})
 *           with r = h − n_-   (normalized homological degree)
 *           bp_from_ranks / bp_add_term → Kh(q,t)
 *
 * GRADING CONVENTION (interim report §3.10–3.11)
 * ----------------------------------------------
 *   r  = normalized homological degree  (cube height h minus n_-)
 *   j  = quantum degree
 *   Kh(q,t) = Σ t^r q^j β_{r,j}
 *   Kh(q,-1) = χ_q = (q + q^{-1}) V_L(q^2)
 *
 *   n_+ = (m + writhe) / 2
 *   n_- = (m − writhe) / 2
 *
 * INTEGER MATRIX CONTRACT
 * -----------------------
 * Topology (or the enhanced-state engine) hands polynomial a Mat that is
 * the boundary d : C_h^{(j)} → C_{h+1}^{(j)} over Z, quantum-preserving.
 *
 * After toSmithForm(m):
 *   free image-rank = # nonzero diagonal entries
 *   torsion coeffs  = |d_ii| whenever |d_ii| > 1
 *
 * Use IntegerMatrix API:
 *   int smith_extract(const Mat *m, int *free_rank,
 *                     int *torsions, int max_tors, int *n_tors);
 *   int toSmithForm_extract(Mat *m, ...);  // mutates m into SNF
 *
 * BIGRADED RANK TABLE CONTRACT
 * ----------------------------
 * Preferred hand-off when topology already reduced the complex:
 *
 *   typedef struct {
 *       int h_min, h_max;     // cube height range (before n_- shift)
 *       int q_min, q_max;     // quantum range
 *       int n_minus;          // for r = h - n_minus
 *       int **chain_dims;     // [h - h_min][q - q_min] = dim C_h^{(q)}
 *       int **image_ranks;    // [h - h_min][q - q_min] = rank(d_h^{(q)})
 *                             // for h in [h_min, h_max); last row unused
 *   } KhRankTable;
 *
 * Polynomial then computes
 *   β(h,q) = chain_dims[h][q] − image_ranks[h][q] − image_ranks[h-1][q]
 * and emits terms at (r, j) = (h - n_minus, q).
 *
 * Helper implemented below:
 *   BivariatePoly *kh_poincare_from_rank_table(const KhRankTable *table);
 *
 * TORSION (optional)
 * ------------------
 * If topology also reports torsion summands on coker(d_h), pass them via
 * KhTorsionReport (see KhPoincare.h). Placement: target degree (h+1, q)
 * with normalized r = (h+1) - n_minus — same as enhanced-state path.
 *
 * WHAT TOPOLOGY MUST GUARANTEE
 * ----------------------------
 * 1. Matrices are over Z with consistent generator order in adjacent degrees.
 * 2. Quantum grading is preserved by each differential (blocks by j).
 * 3. n_+, n_- (or writhe + m) are available for normalization.
 * 4. Deloop / Gauss do not change the homotopy type (homology invariant).
 * 5. For fixed-strand scans, width bounds may be asserted but must not
 *    silently drop generators that affect homology.
 *
 * WHAT POLYNOMIAL GUARANTEES
 * --------------------------
 * 1. kh_poincare / compute_khovanov_polynomial return free Poincaré only.
 * 2. kh_poincare_with_torsion additionally fills KhTorsionReport when asked.
 * 3. bp_eval_t(poly, -1) is the graded Euler characteristic χ_q.
 * 4. Gold tests (unknot, both trefoils) lock the grading convention.
 *
 * NON-GOALS FOR THIS HEADER
 * -------------------------
 * - Implementing FullScanning / deloop (topology/)
 * - Replacing CobMatrix with Mat inside topology
 * - Forcing a single binary; linking is the integrator's job
 */

#include <stdbool.h>
#include "../polynomial/BivariatePoly.h"
#include "KhPoincare.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Bigraded chain dimensions and image ranks before n_- normalization. */
typedef struct {
    int h_min;
    int h_max;       /* inclusive */
    int q_min;
    int q_max;       /* inclusive */
    int n_minus;     /* r = h - n_minus */

    /*
     * chain_dims[hi][qi]  = dim C_{h_min+hi}^{(q_min+qi)}
     * image_ranks[hi][qi] = rank(d_{h_min+hi}^{(q_min+qi)})
     *   defined for hi = 0 .. (h_max - h_min - 1); may be NULL if h_max==h_min
     * Rows must be non-NULL for hi in range where data exists; missing row = 0.
     */
    int **chain_dims;
    int **image_ranks;
} KhRankTable;

/*
 * Assemble free Poincaré polynomial from a rank table.
 * Returns newly allocated BivariatePoly (caller bp_free's).
 * On invalid table returns empty polynomial.
 */
BivariatePoly *kh_poincare_from_rank_table(const KhRankTable *table);

/*
 * Optional: fill torsion report from parallel arrays of the same shape as
 * image_ranks, where tors_counts[hi][qi] and tors_values[hi][qi][k] list
 * |d_ii|>1 on d_h at that block. May be left unimplemented by topology;
 * pass tors=NULL to ignore.
 *
 * Signature reserved for a future topology exporter — not required for
 * the enhanced-state path (which uses kh_poincare_with_torsion directly).
 */
typedef struct {
    int **tors_counts; /* [hi][qi] */
    int ***tors_values; /* [hi][qi][k], k < tors_counts[hi][qi] */
} KhTorsionTable;

BivariatePoly *kh_poincare_from_rank_table_with_torsion(
    const KhRankTable *table,
    const KhTorsionTable *tors_table,
    KhTorsionReport *out_tors);

#ifdef __cplusplus
}
#endif

#endif /* TOPOLOGY_BRIDGE_H */
