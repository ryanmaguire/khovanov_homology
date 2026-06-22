/*
 *  CannedCobordismImpl.h
 *
 *  C translation of org.katlas.JavaKh.CannedCobordismImpl
 *
 *  Purpose:
 *      Defines the concrete implementation data for a CannedCobordism
 *      (a cobordism between two Caps / smoothings). This struct is
 *      stored in CannedCobordism.impl_data and tracks the topological
 *      information: boundary components, connected components, dots,
 *      genus, and the lazy reverse-map arrays.
 */

#ifndef CANNEDCOBORDISMIMPL_H
#define CANNEDCOBORDISMIMPL_H

#include "CannedCobordism.h"
#include "Cap.h"
#include <stdbool.h>
#include <stdint.h>


/*
 *  CannedCobordismImplData
 *
 *  Purpose:
 *      Holds the topological data for a concrete cobordism between
 *      two Caps. This is what gets stored in CannedCobordism.impl_data.
 *
 *  Fields:
 *      n                   — number of boundary points (must match Caps)
 *      hpower              — homological degree offset (power of h)
 *      top, bottom         — pointers to the bounding Caps
 *
 *      nbc                 — number of boundary components
 *      offtop, offbot      — offsets into the connectedComponent array
 *      component           — array[n]: which component each edge belongs to
 *                            (Flat-allocated with the struct)
 *
 *      ncc                 — number of connected components
 *      connectedComponent  — array[nbc]: which connected component each BC belongs to
 *                            (Flat-allocated with the struct)
 *      dots                — array[ncc]: dot count per connected component
 *      genus               — array[ncc]: genus per connected component
 *
 *      (Lazy, computed by reverseMaps:)
 *      boundaryComponents  — array[ncc] of pointers to BC lists
 *      bc_sizes            — array[ncc]: size of each boundaryComponents sub-array
 *      edges               — array[offtop] of pointers to edge lists
 *      edge_sizes          — array[offtop]: size of each edges sub-array
 *      reverse_maps_done   — flag: true after reverseMaps has been called
 *
 *      hashcode            — cached hash code (0 if not yet computed)
 */
typedef struct CannedCobordismImplData {
  int n;
  int hpower;
  Cap *top;
  Cap *bottom;

  int nbc;
  int offtop;
  int offbot;
  int *component;

  int ncc;
  int *connectedComponent;
  int *dots;
  int *genus;

  int *bc_sizes;
  int **boundaryComponents;

  int *edge_sizes;
  int **edges;

  bool reverse_maps_done;

  int hashcode;
} CannedCobordismImplData;

/* ================================================================
 *  Public API
 * ================================================================ */

/*
 *  CannedCobordismImpl_create
 *
 *  Purpose:
 *      Constructs a CannedCobordismImplData by tracing boundary
 *      components through the top and bottom Cap pairings.
 *
 *  Arguments:
 *      top    — the source (top) Cap
 *      bottom — the target (bottom) Cap
 *
 *  Output:
 *      Heap-allocated CannedCobordismImplData. Uses a "Flat Allocation"
 *      strategy: the struct and its primary arrays (component, connectedComponent)
 *      are allocated as a single contiguous memory block to maximize cache
 *      locality and minimize allocation overhead.
 *
 *  Method:
 *      Pre-calculates the maximum possible footprint, allocates once, and
 *      offsets pointers into the block. Then traces edges through pairings
 *      to identify mixed boundary components.
 */
CannedCobordismImplData *CannedCobordismImpl_create(Cap *top, Cap *bottom);

/*
 *  CannedCobordismImpl_isomorphism
 *
 *  Purpose:
 *      Factory: creates the identity/isomorphism cobordism for a Cap.
 *
 *      Ordinary boundary components are represented by cylinder components.
 *      Interior cycle components are represented by components connecting
 *      the corresponding top-cycle and bottom-cycle boundary components.
 *
 *  Arguments:
 *      c — the Cap. May contain interior cycles; each top-cycle boundary is
 *          paired with the corresponding bottom-cycle boundary in the identity
 *          cobordism.
 *
 *  Output:
 *      A fully initialised CannedCobordism (vtable-wrapped) that represents
 *      the identity cobordism on c.
 *
 *  Method:
 *      Build a CannedCobordismImpl with top == bottom == c. Ordinary boundary
 *      components are made into separate cylinder components, and each
 *      interior cycle connects its top and bottom cycle boundaries.
 */
CannedCobordism *CannedCobordismImpl_isomorphism(Cap *c);

/*
 *  CannedCobordismImpl_reverseMaps
 *
 *  Purpose:
 *      Lazily fills the boundaryComponents[][] and edges[][] arrays.
 *      These reverse maps allow looking up which boundary components
 *      and edges belong to each connected component.
 *
 *  Arguments:
 *      impl — the impl data to fill
 *
 *  Output:
 *      None (modifies impl in place).
 *
 *  Method:
 *      Count boundary components per connected component, then
 *      fill the arrays. Count edges per mixed boundary component,
 *      then fill those arrays.
 */
void CannedCobordismImpl_reverseMaps(CannedCobordismImplData *impl);

/*
 *  CannedCobordismImpl_equals
 *
 *  Purpose:
 *      Tests structural equality of two CannedCobordismImpl objects
 *
 *  Arguments:
 *      a — first impl data
 *      b — second impl data
 *
 *  Output:
 *      true if the two impl objects represent the same cobordism data
 *
 *  Method:
 *      Requires matching top and bottom Caps, matching boundary-component
 *      decomposition component[], and matching hpower, dots, and genus.
 *
 *      The connectedComponent[] labels are compared up to relabeling:
 *      two cobordisms are equal if they define the same partition of boundary
 *      components into connected components, even if the internal component
 *      numbers differ
 *
 *      Closed connected components not incident to boundary components are
 *      matched by their dot/genus data
 */
bool CannedCobordismImpl_equals(const CannedCobordismImplData *a,
                                const CannedCobordismImplData *b);

/*
 *  CannedCobordismImpl_hashCode
 *
 *  Purpose:
 *      Computes a hash code matching the Java implementation.
 *
 *  Arguments:
 *      impl — the impl data to hash
 *
 *  Output:
 *      Integer hash code (cached after first computation).
 *
 *  Method:
 *      Combines n, top hash, bottom hash, connectedComponent,
 *      dots, genus, and hpower with bit-shifting.
 */
int CannedCobordismImpl_hashCode(CannedCobordismImplData *impl);

/*
 *  CannedCobordismImpl_compareTo
 *
 *  Purpose:
 *      Total ordering, matching Java's compareTo.
 *
 *  Arguments:
 *      a — first impl data
 *      b — second impl data
 *
 *  Output:
 *      Negative if a < b, 0 if equal, positive if a > b.
 *
 *  Method:
 *      Compare nbc, connectedComponent, ncc, dots, genus, n, hpower,
 *      top/bottom pairings, and ncycles.
 */
int CannedCobordismImpl_compareTo(const CannedCobordismImplData *a,
                                  const CannedCobordismImplData *b);

/*
 *  CannedCobordismImpl_isIsomorphism
 *
 *  Purpose:
 *      Checks whether this cobordism is the identity/isomorphism cobordism
 *      on its Cap.
 *
 *  Arguments:
 *      impl — the impl data to check
 *
 *  Output:
 *      true if the cobordism is the identity/isomorphism cobordism.
 *
 *  Method:
 *      Requires top == bottom, hpower == 0, and no dots or genus.
 *      Each ordinary boundary component must be its own cylinder component.
 *      Each interior cycle must connect exactly the corresponding top-cycle
 *      boundary component to the corresponding bottom-cycle boundary component.
 *
 *      Connected-component labels themselves need not equal boundary indices.
 */
bool CannedCobordismImpl_isIsomorphism(const CannedCobordismImplData *impl);

/*
 *  CannedCobordismImpl_check
 *
 *  Purpose:
 *      Sanity check on the internal consistency of the impl data.
 *
 *  Arguments:
 *      impl — the impl data to check
 *
 *  Output:
 *      true if consistent.
 *
 *  Method:
 *      Verify nbc/ncc bounds and connectedComponent values.
 */
bool CannedCobordismImpl_check(const CannedCobordismImplData *impl);

/*
 *  CannedCobordismImpl_compose_vertical
 *
 *  Purpose:
 *      Vertical composition: stacks cc on top of this. The top of
 *      this must equal the bottom of cc.
 *
 *  Arguments:
 *      this_impl — bottom cobordism
 *      cc_impl   — top cobordism (cc goes on top of this)
 *
 *  Output:
 *      A new CannedCobordism (vtable-wrapped) representing the
 *      composition.
 *
 *  Method:
 *      1. Call reverseMaps on both inputs.
 *      2. Create a new impl via Flat Allocation.
 *      3. Merge connected components using stack-allocated VLAs (safe due to
 *         Fixed-Strand Conjecture bounds).
 *      4. Compute genus via Euler characteristic.
 *      5. Relabel connected components in ascending order.
 */
CannedCobordism *
CannedCobordismImpl_compose_vertical(CannedCobordismImplData *this_impl,
                                     CannedCobordismImplData *cc_impl);

/*
 *  CannedCobordismImpl_compose_horizontal
 *
 *  Purpose:
 *      Horizontal (partial) composition: joins nc consecutive edges
 *      between this and cc.
 *
 *  Arguments:
 *      this_impl — first cobordism
 *      start     — starting edge index in this
 *      cc_impl   — second cobordism
 *      cstart    — starting edge index in cc
 *      nc        — number of edges to join
 *
 *  Output:
 *      A new CannedCobordism (vtable-wrapped) representing the
 *      horizontal composition.
 *
 *  Method:
 *      1. Compose the top and bottom Caps via Cap_compose.
 *      2. Create new impl from the composed caps.
 *      3. Merge connected components, tracking middle join edges
 *         with midConComp[].
 *      4. Handle newly created cycles from tjoins/bjoins.
 *      5. Compute genus via Euler characteristic.
 *      6. Relabel and sort connected components.
 */
CannedCobordism *CannedCobordismImpl_compose_horizontal(
    CannedCobordismImplData *this_impl, int start,
    CannedCobordismImplData *cc_impl, int cstart, int nc);

/*
 *  CannedCobordismImpl_as_CannedCobordism
 *
 *  Purpose:
 *      Wraps a CannedCobordismImplData in a CannedCobordism vtable
 *      struct so it can be used polymorphically.
 *
 *  Arguments:
 *      impl — the impl data to wrap
 *
 *  Output:
 *      A heap-allocated CannedCobordism with vtable function pointers
 *      pointing to CannedCobordismImpl dispatch functions.
 *
 *  Method:
 *      Allocate a CannedCobordism, set source/target from impl,
 *      store impl in impl_data, wire vtable pointers.
 */
CannedCobordism *
CannedCobordismImpl_as_CannedCobordism(CannedCobordismImplData *impl);

/*
 *  CannedCobordismImpl_free
 *
 *  Purpose:
 *      Frees all memory owned by a CannedCobordismImplData.
 *
 *  Arguments:
 *      impl — the impl data to free (may be NULL)
 *
 *  Output:
 *      None.
 *
 *  Method:
 *      Frees the primary Flat Allocation block (which includes component and
 *      connectedComponent) and then frees the lazily-allocated reverse mapping
 *      structures.
 */
void CannedCobordismImpl_free(CannedCobordismImplData *impl);
/*
 *  CannedCobordismImpl_stripClosed
 *
 *  Purpose:
 *      Remove closed connected components that have already been accounted
 *      for by reduction/evaluation, so that remaining cobordisms compare by
 *      their boundary-attached data.
 *
 *  Warning:
 *      Do not call this before applying the Bar-Natan closed-surface rules;
 *      closed components with dots or genus may contribute nontrivial scalar
 *      factors.
 */
CannedCobordism *CannedCobordismImpl_stripClosed(CannedCobordism *cc);
CannedCobordism *CannedCobordismImpl_capOffTop(CannedCobordism *cc, Cap *new_top, bool add_dot);
CannedCobordism *CannedCobordismImpl_cupOnBottom(CannedCobordism *cc, Cap *new_bottom, bool add_dot);
#endif /* CANNEDCOBORDISMIMPL_H */
