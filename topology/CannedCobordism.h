/*
 * Purpose:
 *      Represents a canned cobordism between Caps. By using a struct containing
 *      function pointers (vtable), concrete implementations can provide their
 * own behaviour.
 */

#ifndef CANNEDCOBORDISM_H
#define CANNEDCOBORDISM_H

#include <stdbool.h>
#include <stddef.h>

#include "Cap.h"
typedef struct CannedCobordism CannedCobordism;

/*
 * Purpose:
 *      Abstract interface struct for a CannedCobordism.
 * Arguments:
 *      None
 * Argument descriptions:
 *      None
 * Output:
 *      None
 * Output description:
 *      None
 * Method:
 *      Contains fields for source and target caps, and a vtable of function
 * pointers for composition, partial composition, isomorphism checks, map
 * reversal, and memory cleanup.
 */
struct CannedCobordism {
  Cap *source;
  Cap *target;
  void *impl_data;

  /* vtable function: full composition */
  CannedCobordism *(*compose)(const CannedCobordism *self,
                              const CannedCobordism *other);

  /* vtable function: partial composition */
  CannedCobordism *(*compose_partial)(const CannedCobordism *self, int start,
                                      const CannedCobordism *cc, int cstart,
                                      int nc);

  /* vtable function: check isomorphism */
  bool (*isIsomorphism)(const CannedCobordism *self);

  /* vtable function: reverse maps */
  void (*reverseMaps)(CannedCobordism *self);

  /* vtable function: free implementation */
  void (*free_impl)(CannedCobordism *self);
};

/*
 * Purpose:
 *      These dispatch through the function pointers so callers don't
 *      need to know the vtable layout.
 */

static inline Cap *CannedCobordism_source(const CannedCobordism *cc) {
  return cc->source;
}

static inline Cap *CannedCobordism_target(const CannedCobordism *cc) {
  return cc->target;
}

static inline CannedCobordism *
CannedCobordism_compose(const CannedCobordism *cc,
                        const CannedCobordism *other) {
  return cc->compose(cc, other);
}

static inline CannedCobordism *
CannedCobordism_compose_partial(const CannedCobordism *cc, int start,
                                const CannedCobordism *other, int cstart,
                                int nc) {
  return cc->compose_partial(cc, start, other, cstart, nc);
}

static inline bool CannedCobordism_isIsomorphism(const CannedCobordism *cc) {
  return cc->isIsomorphism(cc);
}

static inline void CannedCobordism_reverseMaps(CannedCobordism *cc) {
  cc->reverseMaps(cc);
}

static inline void CannedCobordism_free(CannedCobordism *cc) {
  if (cc == NULL)
    return;
  if (cc->free_impl)
    cc->free_impl(cc);
}

#endif /* CANNEDCOBORDISM_H */
