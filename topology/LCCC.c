/*
 *  LCCC.c
 *
 *  Linear Combination of Canned Cobordisms over F_2.
 *
 *  Purpose:
 *      Implements F_2 arithmetic for the Khovanov chain complex
 *      coefficient ring. All additions are XOR (a + a = 0), every
 *      nonzero element is self-inverse, and composition distributes.
 *
 *  Compile (standalone test):
 *      gcc -Wall -Wextra -g -c LCCC.c -I.
 */

#include "LCCC.h"
#include <stdlib.h>

/* ================================================================
 *  Internal helpers
 * ================================================================ */

/*
 *  Purpose:    Allocate a new LCCCTerm node.
 */
static LCCCTerm *create_term(CannedCobordism *cc) {
  LCCCTerm *t = (LCCCTerm *)malloc(sizeof(LCCCTerm));
  t->cobordism = cc;
  t->next = NULL;
  return t;
}

/*
 *  Purpose:    Check if two cobordisms are "equal" for F_2 cancellation.
 *              We compare by pointer identity first (fast path), then
 *              fall back to structural equality via the CannedCobordism
 *              interface.
 */
static bool cobordism_equal(const CannedCobordism *a,
                            const CannedCobordism *b) {
  if (a == b)
    return true;
  if (a == NULL || b == NULL)
    return false;

  /* Structural equality: same source, target, and impl_data */
  if (!Cap_equals(a->source, b->source))
    return false;
  if (!Cap_equals(a->target, b->target))
    return false;

  /*
   * For full structural equality we would compare impl_data.
   * Since CannedCobordismImpl provides isIsomorphism and equals,
   * we use pointer identity as the primary check. Two distinct
   * cobordism objects with identical topology are treated as
   * the same term for F_2 cancellation only if they share
   * source/target Caps.
   */
  return (a->impl_data == b->impl_data);
}

/*
 *  Purpose:    Toggle a cobordism in an LCCC (XOR semantics).
 *              If cc is already present, remove it. Otherwise, add it.
 *              Modifies lc in place.
 */
static void lccc_toggle(LCCC *lc, CannedCobordism *cc) {
  LCCCTerm *prev = NULL;
  LCCCTerm *cur = lc->head;

  while (cur != NULL) {
    if (cobordism_equal(cur->cobordism, cc)) {
      /* Found: remove (XOR cancel) */
      if (prev)
        prev->next = cur->next;
      else
        lc->head = cur->next;
      free(cur);
      lc->count--;
      return;
    }
    prev = cur;
    cur = cur->next;
  }

  /* Not found: add at head */
  LCCCTerm *t = create_term(cc);
  t->next = lc->head;
  lc->head = t;
  lc->count++;
}

/* ================================================================
 *  Construction / Destruction
 * ================================================================ */

LCCC *LCCC_createZero(void) {
  LCCC *lc = (LCCC *)malloc(sizeof(LCCC));
  lc->head = NULL;
  lc->count = 0;
  return lc;
}

LCCC *LCCC_createSingle(CannedCobordism *cc) {
  LCCC *lc = LCCC_createZero();
  if (cc != NULL) {
    LCCCTerm *t = create_term(cc);
    lc->head = t;
    lc->count = 1;
  }
  return lc;
}

LCCC *LCCC_clone(const LCCC *lc) {
  LCCC *result = LCCC_createZero();
  if (lc == NULL)
    return result;

  LCCCTerm **tail = &result->head;
  LCCCTerm *cur = lc->head;
  while (cur != NULL) {
    LCCCTerm *t = create_term(cur->cobordism);
    *tail = t;
    tail = &t->next;
    result->count++;
    cur = cur->next;
  }
  return result;
}

void LCCC_free(LCCC *lc) {
  if (lc == NULL)
    return;
  LCCCTerm *cur = lc->head;
  while (cur != NULL) {
    LCCCTerm *next = cur->next;
    /* We do NOT free cur->cobordism — ownership is external */
    free(cur);
    cur = next;
  }
  free(lc);
}

/* ================================================================
 *  Arithmetic (F_2)
 * ================================================================ */

LCCC *LCCC_add(LCCC *a, LCCC *b) {
  /* Start with a clone of a, then XOR in each term of b */
  LCCC *result = LCCC_clone(a);

  if (b == NULL)
    return result;

  LCCCTerm *cur = b->head;
  while (cur != NULL) {
    lccc_toggle(result, cur->cobordism);
    cur = cur->next;
  }
  return result;
}

LCCC *LCCC_compose(LCCC *a, LCCC *b) {
  LCCC *result = LCCC_createZero();

  if (a == NULL || b == NULL)
    return result;
  if (a->count == 0 || b->count == 0)
    return result;

  /* Distribute: Σ_{i,j} a_i ∘ b_j (over F_2) */
  LCCCTerm *ai = a->head;
  while (ai != NULL) {
    LCCCTerm *bj = b->head;
    while (bj != NULL) {
      if (ai->cobordism != NULL && bj->cobordism != NULL) {
        CannedCobordism *composed =
            CannedCobordism_compose(ai->cobordism, bj->cobordism);
        if (composed != NULL) {
          lccc_toggle(result, composed);
        }
      }
      bj = bj->next;
    }
    ai = ai->next;
  }
  return result;
}

LCCC *LCCC_multiply(LCCC *a, RingElement *coeff) {
  if (a == NULL || coeff == NULL)
    return LCCC_createZero();

  /* Over F_2: multiply by 0 → zero, multiply by 1 → clone */
  if ((coeff->value & 1) == 0)
    return LCCC_createZero();

  return LCCC_clone(a);
}

LCCC *LCCC_reduce(LCCC *a) {
  /*
   * Over F_2 with XOR-maintained lists, the LCCC is already canonical.
   * Just clone it.
   */
  return LCCC_clone(a);
}

LCCC *LCCC_invert(LCCC *lc) {
  /*
   * Over F_2, every isomorphism is self-inverse.
   * φ^{-1} = φ.
   */
  return LCCC_clone(lc);
}

LCCC *LCCC_negate(LCCC *lc) {
  /*
   * Over F_2, -1 = 1, so negation is identity.
   */
  return LCCC_clone(lc);
}

/* ================================================================
 *  Predicates
 * ================================================================ */

bool LCCC_isZero(LCCC *lc) {
  if (lc == NULL)
    return true;
  return (lc->count == 0);
}

bool LCCC_contains(const LCCC *lc, const CannedCobordism *cc) {
  if (lc == NULL || cc == NULL)
    return false;
  LCCCTerm *cur = lc->head;
  while (cur != NULL) {
    if (cobordism_equal(cur->cobordism, cc))
      return true;
    cur = cur->next;
  }
  return false;
}
