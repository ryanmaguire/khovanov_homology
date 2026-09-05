#ifndef KNOT_ENCODER_H
#define KNOT_ENCODER_H

#include "../polynomial/BivariatePoly.h"

/* knot representation: planar diagram with 0- or 1-based labels depending
 * on construction path. knot_from_pd() converts 1-based input to 0-based.
 * Enhanced-state homology (kh_poincare) expects 0-based labels in 0..2m-1. */
typedef struct {
    int **pd;   /* planar diagram codes: m crossings × 4 edge labels */
    int writhe; /* writhe of the oriented diagram */
    int *edges; /* optional 4m edge labels (may be NULL) */
    int m;      /* number of crossings */
} Knot;

/* resolution register |r> (m bits) — reserved for future use */
typedef struct {
    int *r;
    int m;
} Resolution;

/* resolved knot register |K_r> */
typedef struct {
    int *resolved_edges;
    int m;
} ResolvedKnot;

/* loop enumeration register |L_r> */
typedef struct {
    int *loop_min_labels;
    int num_loops;
} LoopRegister;

/* enhanced state |s> (labels on loops: 0=1, 1=X) */
typedef struct {
    int *s;
    int num_loops;
} EnhancedState;

/* hard-coded positive trefoil via knot_from_pd (1-based PD → 0-based) */
Knot *knot_create_trefoil(void);

/* create knot from 1-based PD code array; stores 0-based labels internally */
Knot *knot_from_pd(int pd[][4], int m);

/* count circles in resolution r (bitmask), using union-find on edges */
int knot_count_loops(Knot *knot, int r);

void knot_free(Knot *k);

/*
 * TRUE Khovanov / Poincaré polynomial Kh(q,t) = Σ t^r q^j β_{r,j}.
 *
 * Implementation: delegates to kh_poincare() (enhanced states + Smith form).
 * This is the production entry point for free ranks.
 *
 * Previously this name incorrectly returned a signed state sum (graded Euler
 * characteristic of the cube, not homology). That behaviour is now exposed
 * only as compute_graded_euler_state_sum().
 */
BivariatePoly *compute_khovanov_polynomial(Knot *knot);

/*
 * Graded Euler / signed state-sum over resolutions (NOT homology).
 *
 * For each resolution r with h = popcount(r) and ℓ loops:
 *   add  (-1)^h  at  q^{w + h - 2ℓ}  t^h
 * (cube grading; not normalized by n_±).
 *
 * Kept for debugging / comparison with the Kauffman bracket path.
 * Do not use this as a substitute for Kh(q,t).
 */
BivariatePoly *compute_graded_euler_state_sum(Knot *knot);

#endif /* KNOT_ENCODER_H */
