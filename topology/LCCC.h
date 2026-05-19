/*
 *  LCCC.h
 *
 *  Linear Combination of Canned Cobordisms over F_2.
 *
 *  Purpose:
 *      Provides the concrete algebraic coefficient type for the
 *      Khovanov chain complex. Over F_2 (mod 2), coefficients are
 *      implicit: a term is either present (coefficient 1) or absent
 *      (coefficient 0). Addition is XOR, and every nonzero morphism
 *      is its own inverse.
 *
 *  Structure:
 *      An LCCC is a singly-linked list of CannedCobordism pointers.
 *      Each node represents a term with implicit coefficient 1.
 *      Duplicate terms cancel (XOR semantics).
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
  CannedCobordism *cobordism;
  struct LCCCTerm *next;
} LCCCTerm;

/*
 *  LCCC — the linear combination itself.
 */
typedef struct LCCC {
  LCCCTerm *head;
  int count; /* Number of terms (for quick isZero check) */
} LCCC;

/*
 *  RingElement — coefficient ring element.
 *  Over F_2 this is trivial (0 or 1), but we keep the type
 *  for API compatibility with CobMatrix.h forward declarations.
 */
typedef struct RingElement {
  int value; /* 0 or 1 over F_2 */
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
 *  Output:     Heap-allocated LCCC with one term.
 */
LCCC *LCCC_createSingle(CannedCobordism *cc);

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
 *  Arithmetic (F_2)
 * ================================================================ */

/*
 *  Purpose:    Add two LCCCs over F_2 (XOR semantics).
 *              Creates a new LCCC; inputs are not modified.
 *  Method:     Merge term lists. If a cobordism appears in both,
 *              it cancels (removed from result).
 */
LCCC *LCCC_add(LCCC *a, LCCC *b);

/*
 *  Purpose:    Compose two LCCCs: distribute composition over all
 *              pairs of terms and sum (over F_2).
 *              (a1 + a2 + ...) ∘ (b1 + b2 + ...) =
 *              Σ_{i,j} ai ∘ bj
 */
LCCC *LCCC_compose(LCCC *a, LCCC *b);

/*
 *  Purpose:    Multiply an LCCC by a ring element.
 *              Over F_2: multiply by 0 gives zero, by 1 gives clone.
 */
LCCC *LCCC_multiply(LCCC *a, RingElement *coeff);

/*
 *  Purpose:    Reduce an LCCC (combine like terms).
 *              Over F_2, this is a no-op on a well-formed LCCC
 *              but we canonicalize just in case.
 *  Output:     New reduced LCCC (caller must free).
 */
LCCC *LCCC_reduce(LCCC *a);

/*
 *  Purpose:    Invert a unit LCCC (the pivot isomorphism φ).
 *              Over F_2, every nonzero LCCC that is an isomorphism
 *              is its own inverse: φ^{-1} = φ.
 *  Output:     Clone of lc.
 */
LCCC *LCCC_invert(LCCC *lc);

/*
 *  Purpose:    Negate an LCCC.
 *              Over F_2, negation is identity.
 *  Output:     Clone of lc.
 */
LCCC *LCCC_negate(LCCC *lc);

/* ================================================================
 *  Predicates
 * ================================================================ */

/*
 *  Purpose:    Check if the LCCC is the zero element (empty list).
 */
bool LCCC_isZero(LCCC *lc);

/*
 *  Purpose:    Check if a single cobordism term is present in the LCCC.
 */
bool LCCC_contains(const LCCC *lc, const CannedCobordism *cc);

#endif /* LCCC_H */
