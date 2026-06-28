/*
 *  LCCC.h
 *
 *  Linear Combination of Canned Cobordisms with integer coefficients.
 *
 *  Purpose:
 *      Provides the concrete algebraic coefficient type for the
 *      Khovanov chain complex. Each term stores an explicit integer
 *      coefficient together with a canned cobordism.
 *
 *  Structure:
 *      An LCCC is a singly-linked list of (coefficient, cobordism)
 *      pairs. Structurally equal cobordisms are combined by adding
 *      coefficients, and zero coefficients are removed.
 */

#ifndef LCCC_H
#define LCCC_H

#include "CannedCobordism.h"
#include <stdbool.h>

/* ================================================================
 *  Type Definitions
 * ================================================================ */

/*
 *  LCCCTerm — a single node in the linked list.
 */
typedef struct LCCCTerm {
  int coeff; /* Integer coefficient of this cobordism term. */
  CannedCobordism *cobordism;
  struct LCCCTerm *next;
} LCCCTerm;

/*
 *  LCCC — the linear combination itself.
 */
typedef struct LCCC {
  LCCCTerm *head;
  int count; /* Number of terms, for quick zero checks. */
} LCCC;

/*
 *  RingElement — coefficient ring element.
 *  We keep this lightweight wrapper for API compatibility with
 *  CobMatrix.h forward declarations.
 */
typedef struct RingElement {
  int value;
} RingElement;

/* ================================================================
 *  Construction / Destruction
 * ================================================================ */

/*
 *  Purpose:    Create an empty LCCC (the zero element).
 *  Output:     Heap-allocated LCCC with no terms.
 */
LCCC *LCCC_createZero(void);

/*
 *  Purpose:    Create an LCCC containing a single cobordism term.
 *  Arguments:  cc — the cobordism (ownership is NOT transferred;
 *              the LCCC stores a pointer to cc).
 *              coeff — integer coefficient for the term.
 *  Output:     Heap-allocated LCCC with one term when coeff != 0.
 */
LCCC *LCCC_createSingle(CannedCobordism *cc, int coeff);

/*
 *  Purpose:    Deep-clone an LCCC (copies the list structure,
 *              but cobordism pointers are shared).
 */
LCCC *LCCC_clone(const LCCC *lc);

/*
 *  Purpose:    Free an LCCC and its term list.
 *              Does NOT free the cobordism pointers themselves.
 */
void LCCC_free(LCCC *lc);

/* ================================================================
 *  Arithmetic
 * ================================================================ */

/*
 *  Purpose:    Add two LCCCs by combining coefficients of equal terms.
 *              Creates a new LCCC; inputs are not modified.
 *  Method:     Merge term lists and remove any term whose resulting
 *              coefficient becomes zero.
 */
LCCC *LCCC_add(LCCC *a, LCCC *b);

/*
 *  Purpose:    Compose two LCCCs: distribute composition over all
 *              pairs of terms and sum with integer coefficients.
 *              (a1 + a2 + ...) ∘ (b1 + b2 + ...) =
 *              Σ_{i,j} ai ∘ bj
 */
LCCC *LCCC_compose(LCCC *a, LCCC *b);

/*
 *  Purpose:    Multiply an LCCC by an integer ring element.
 */
LCCC *LCCC_multiply(LCCC *a, RingElement *coeff);

/*
 *  Purpose:    Reduce an LCCC using the Bar-Natan surface relations.
 *              Closed components are evaluated using the sphere,
 *              dotted-sphere, double-dot, and torus relations.
 *              Components with multiple boundary components are
 *              expanded using neck-cutting.
 *  Output:     New reduced LCCC (caller must free).
 */
LCCC *LCCC_reduce(LCCC *a);

/*
 *  Purpose:    Invert a unit LCCC used as a Gaussian-elimination pivot.
 *              Only single-term pivots with coefficient +/-1 are
 *              invertible over Z. The current elimination path uses
 *              identity-like cobordism pivots, whose inverse is
 *              represented by the same underlying cobordism term.
 *  Output:     Clone of lc when lc is an invertible unit pivot;
 *              NULL otherwise.
 */
LCCC *LCCC_invert(LCCC *lc);

/*
 *  Purpose:    Negate an LCCC by multiplying every coefficient by -1.
 *  Output:     New negated LCCC.
 */
LCCC *LCCC_negate(LCCC *lc);

LCCC *LCCC_capOffTop(LCCC *lc, Cap *new_top, bool add_dot);
LCCC *LCCC_cupOnBottom(LCCC *lc, Cap *new_bottom, bool add_dot);


/*
 *  Purpose:    Check if the LCCC is the zero element (empty list).
 */
bool LCCC_isZero(LCCC *lc);

/*
 *  Purpose:    Check if a single cobordism term is present in the LCCC.
 */
bool LCCC_contains(const LCCC *lc, const CannedCobordism *cc);

#endif //LCCC_H
