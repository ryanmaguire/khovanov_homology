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
 *      n                   — number of boundary edges (same as top->n)
 *      hpower              — power of h (delooping parameter)
 *      top                 — source Cap (top boundary)
 *      bottom              — target Cap (bottom boundary)
 *      nbc                 — total number of boundary components
 *                            (mixed + top cycles + bottom cycles)
 *      offtop              — offset: boundary components [0, offtop) are
 *                            mixed; [offtop, offbot) are top cycles
 *      offbot              — offset: [offbot, nbc) are bottom cycles
 *      component           — array[n]: which mixed boundary component
 *                            each edge belongs to
 *      ncc                 — number of connected components
 *      connectedComponent  — array[nbc]: which connected component
 *                            each boundary component belongs to
 *      dots                — array[ncc]: dot count per connected component
 *      genus               — array[ncc]: genus per connected component
 *
 *      (Lazy, computed by reverseMaps:)
 *      boundaryComponents  — array[ncc] of arrays: which boundary
 *                            components belong to each connected component
 *      bc_sizes            — array[ncc]: size of each boundaryComponents
 *                            sub-array
 *      edges               — array[offtop] of arrays: which edges belong
 *                            to each mixed boundary component
 *      edge_sizes          — array[offtop]: size of each edges sub-array
 *      reverse_maps_done   — flag: true after reverseMaps has been called
 *      hashcode            — cached hash code (0 if not yet computed)
 */
typedef struct CannedCobordismImplData {
  int n;
  int hpower;
  Cap *top;
  Cap *bottom;

  int8_t nbc;
  int offtop;
  int offbot;
  int8_t *component; /* array[n]   */

  int8_t ncc;
  int8_t *connectedComponent; /* array[nbc] */
  int8_t *dots;               /* array[ncc] */
  int8_t *genus;              /* array[ncc] */

  /* Lazy reverse maps (filled by CannedCobordismImpl_reverseMaps) */
  int8_t **boundaryComponents; /* array[ncc] of int8_t arrays */
  int *bc_sizes;               /* array[ncc] */
  int8_t **edges;              /* array[offtop] of int8_t arrays */
  int *edge_sizes;             /* array[offtop] */
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
 *      Heap-allocated CannedCobordismImplData with component[]
 *      and connectedComponent[] initialised. The caller must fill
 *      in ncc, connectedComponent mappings, dots, and genus
 *      (or use a factory function that does so).
 *
 *  Method:
 *      Trace edges through top.pairings then bottom.pairings to
 *      identify mixed boundary components. Set offtop/offbot
 *      from the Cap ncycles values.
 */
CannedCobordismImplData *CannedCobordismImpl_create(Cap *top, Cap *bottom);

/*
 *  CannedCobordismImpl_isomorphism
 *
 *  Purpose:
 *      Factory: creates the identity/isomorphism cobordism for a Cap.
 *      The resulting cobordism has ncc == nbc, each boundary component
 *      is its own connected component, with zero dots and genus.
 *
 *  Arguments:
 *      c — the Cap (must have ncycles == 0)
 *
 *  Output:
 *      A fully initialised CannedCobordism (vtable-wrapped) that
 *      is an isomorphism.
 *
 *  Method:
 *      Build a CannedCobordismImpl with top == bottom == c,
 *      connectedComponent[i] = i, dots = genus = zeros.
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
 *      Tests structural equality of two CannedCobordismImpl objects.
 *
 *  Arguments:
 *      a — first impl data
 *      b — second impl data
 *
 *  Output:
 *      true if all fields match.
 *
 *  Method:
 *      Compare ncc, dots, genus, nbc, connectedComponent, n, hpower,
 *      component, top, and bottom.
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
 *      Checks whether this cobordism is an isomorphism.
 *      Only works on delooped cobordisms (top.ncycles == 0).
 *
 *  Arguments:
 *      impl — the impl data to check
 *
 *  Output:
 *      true if the cobordism is the identity.
 *
 *  Method:
 *      Verify top == bottom, ncycles == 0, nbc == ncc, hpower == 0,
 *      connectedComponent[i] == i, dots[i] == 0, genus[i] == 0.
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
 *      2. Create a new impl with top = cc.top, bottom = this.bottom.
 *      3. Merge connected components from both, tracking middle
 *         cycles with midConComp[].
 *      4. Compute genus via Euler characteristic.
 *      5. Relabel connected components in ascending order.
 *      6. Collect unconnected components (genus/dots).
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
 *      Free all owned arrays (component, connectedComponent, dots,
 *      genus, boundaryComponents, edges, bc_sizes, edge_sizes),
 *      then free the struct itself.
 *      Note: top and bottom Caps are NOT freed (not owned).
 */
void CannedCobordismImpl_free(CannedCobordismImplData *impl);

#endif /* CANNEDCOBORDISMIMPL_H */
