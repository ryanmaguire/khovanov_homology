/*
 * Purpose:
 *      A concrete CannedCobordism that acts as the identity morphism:
 *        - source() and target() return the same Cap
 *        - compose(other)  returns 'other' unchanged
 *        - compose_partial delegates to the other cobordism
 *        - isIsomorphism() always returns true
 *        - reverseMaps()   is a no-op
 */

#include "IdentityCannedCobordism.h"
#include <stdlib.h>

/* ----------------------------------------------------------------
 *  Vtable implementations (static – internal linkage)
 * ---------------------------------------------------------------- */

/*
 * Purpose:
 *      Identity compose: composing with the identity just returns the other.
 * Arguments:
 *      self, other
 * Argument descriptions:
 *      self  - The identity cobordism
 *      other - The other cobordism to compose with
 * Output:
 *      CannedCobordism pointer
 * Output description:
 *      Returns the 'other' cobordism unchanged.
 * Method:
 *      Simply return the 'other' pointer. It represents a shallow "return"
 * without cloning.
 */
static CannedCobordism *identity_compose(const CannedCobordism *self,
                                         const CannedCobordism *other) {
  (void)self; /* unused */
  /* Return the other cobordism unchanged.
   * Note: this is a shallow "return" – we are not cloning.
   * Matches the Java semantics exactly. */
  return (CannedCobordism *)other;
}

/*
 * Purpose:
 *      Partial compose: delegates to the other cobordism with swapped roles.
 * Arguments:
 *      self, start, cc, cstart, nc
 * Argument descriptions:
 *      self   - The identity cobordism
 *      start  - offset in identity cobordism
 *      cc     - The other cobordism
 *      cstart - offset in the other cobordism
 *      nc     - number of edges to join
 * Output:
 *      CannedCobordism pointer
 * Output description:
 *      Returns the result of the composition.
 * Method:
 *      Delegates to cc->compose_partial(cc, cstart, self, start, nc).
 */
static CannedCobordism *identity_compose_partial(const CannedCobordism *self,
                                                 int start,
                                                 const CannedCobordism *cc,
                                                 int cstart, int nc) {
  /* Delegate: cc.compose_partial(cstart, self, start, nc) */
  return cc->compose_partial(cc, cstart, self, start, nc);
}

/*
 * Purpose:
 *      Indicates whether this cobordism is an isomorphism.
 * Arguments:
 *      self
 * Argument descriptions:
 *      self - The identity cobordism
 * Output:
 *      boolean
 * Output description:
 *      Always an isomorphism, so returns true.
 * Method:
 *      Returns true.
 */
static bool identity_isIsomorphism(const CannedCobordism *self) {
  (void)self;
  return true;
}

/*
 * Purpose:
 *      Reverses the internal maps of this cobordism.
 * Arguments:
 *      self
 * Argument descriptions:
 *      self - The identity cobordism
 * Output:
 *      None
 * Output description:
 *      None
 * Method:
 *      No-op for the identity cobordism.
 */
static void identity_reverseMaps(CannedCobordism *self) { (void)self; }

/*
 * Purpose:
 *      Frees implementation-specific data.
 * Arguments:
 *      self
 * Argument descriptions:
 *      self - The identity cobordism
 * Output:
 *      None
 * Output description:
 *      None
 * Method:
 *      The identity cobordism has no extra impl_data to clean up.
 *      We just free the struct itself. Note that the Cap is NOT owned
 *      by the identity cobordism.
 */
static void identity_free_impl(CannedCobordism *self) { free(self); }

/*
 * Purpose:
 *      Constructor for the identity cobordism.
 * Arguments:
 *      cap
 * Argument descriptions:
 *      cap - The Cap to associate with this identity cobordism.
 * Output:
 *      CannedCobordism pointer
 * Output description:
 *      The newly allocated identity cobordism.
 * Method:
 *      Allocate memory, set the capability fields, and wire up the function
 * pointers.
 */

CannedCobordism *IdentityCannedCobordism_create(Cap *cap) {
  CannedCobordism *cc = (CannedCobordism *)malloc(sizeof(CannedCobordism));
  if (cc == NULL)
    return NULL;

  /* Source and target are the same cap (identity morphism) */
  cc->source = cap;
  cc->target = cap;
  cc->impl_data = NULL; /* No extra data needed */

  /* Wire up the vtable */
  cc->compose = identity_compose;
  cc->compose_partial = identity_compose_partial;
  cc->isIsomorphism = identity_isIsomorphism;
  cc->reverseMaps = identity_reverseMaps;
  cc->free_impl = identity_free_impl;

  return cc;
}
