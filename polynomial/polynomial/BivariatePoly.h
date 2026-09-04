#ifndef BIVARIATE_POLY_H
#define BIVARIATE_POLY_H
#include <stdbool.h>

/*
 * Bivariate polynomial for the Khovanov / Poincaré polynomial
 *
 *   Kh(L)(q, t) = sum_{r,j}  t^r q^j  dim Kh^{r,j}(L)
 *
 * Convention (interim report §3.10–3.11):
 *   Setting t = -1 recovers the graded Euler characteristic
 *     χ_q(Kh(L)) = (q + q^{-1}) V_L(q^2)
 *   where V_L is the Jones polynomial in the usual q-variable after
 *   the standard substitution.
 *
 * Monomial exponents:
 *   q_exp  = quantum (j) degree
 *   t_exp  = homological (r) degree  (already normalized, e.g. h - n_-)
 *   coeff  = free rank (integer; torsion is tracked separately upstream)
 */

typedef struct {
    int q_exp; /* quantum degree j */
    int t_exp; /* homological degree r */
    int coeff;
} Monomial;

typedef struct {
    Monomial *terms; /* sparse array of monomials */
    int num_terms;   /* current number of terms */
    int capacity;    /* allocated size */
} BivariatePoly;

/* ---- lifecycle ---- */
BivariatePoly *bp_create(void);
void bp_free(BivariatePoly *p);
BivariatePoly *bp_copy(const BivariatePoly *p);

/* ---- status ---- */
bool bp_is_zero(const BivariatePoly *p);
bool bp_equals(const BivariatePoly *p1, const BivariatePoly *p2);

/* ---- mutation ---- */
/* Merge like terms; zero coefficients are dropped. Keeps terms sorted
 * by (t_exp ascending, then q_exp ascending). */
void bp_add_term(BivariatePoly *p, int q, int t, int c);

/* In-place shift: every term (q,t) -> (q+dq, t+dt). Used for n_±
 * normalization when ranks arrive in unnormalized cube grading. */
void bp_shift(BivariatePoly *p, int dq, int dt);

/* Force canonical order (t then q). Usually unnecessary if only
 * bp_add_term is used, but safe to call before print/compare. */
void bp_sort(BivariatePoly *p);

/* ---- algebra (allocate new results; caller must free) ---- */
BivariatePoly *bp_add(const BivariatePoly *p1, const BivariatePoly *p2);
BivariatePoly *bp_multiply(const BivariatePoly *p1, const BivariatePoly *p2);

/* Evaluate at a fixed integer value of t. Returns a new polynomial
 * in q only (all t_exp == 0). Primary use: t0 = -1 recovers χ_q. */
BivariatePoly *bp_eval_t(const BivariatePoly *p, int t0);

/* ---- assembly from homology ranks ---- */
/*
 * Build Kh(q,t) from a 2D table of free Betti numbers.
 *
 *   betti[h_index][q_index] = dim Kh^{r,j}
 *   h runs h_min .. h_max  (inclusive), already the true homological r
 *   q runs q_min .. q_max  (inclusive), already the true quantum j
 *
 * Rows outside [0, h_max-h_min] or null pointers are treated as zero.
 * Only strictly positive ranks become terms.
 */
BivariatePoly *bp_from_ranks(int **betti,
                             int h_min, int h_max,
                             int q_min, int q_max);

/* ---- output ---- */
void bp_print(const BivariatePoly *p);
/* Heap-allocated C string; caller must free(). */
char *bp_to_string(const BivariatePoly *p);

#endif /* BIVARIATE_POLY_H */
