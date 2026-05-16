/*
 *  CannedCobordismImpl.c
 *
 *
 *  Purpose:
 *      Concrete implementation of CannedCobordism representing a
 *      cobordism surface between two Caps in the Khovanov complex.
 */

#include "CannedCobordismImpl.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 *  Static utility arrays (matching Java's zeros/counting)
 * ================================================================ */

#define ZEROSIZE 127

static int8_t zero_buf[ZEROSIZE]; /* all zeros, shared */

/* ================================================================
 *  Vtable dispatch functions (static)
 * ================================================================ */

static CannedCobordism *impl_compose(const CannedCobordism *self,
                                     const CannedCobordism *other) {
  CannedCobordismImplData *sd = (CannedCobordismImplData *)self->impl_data;
  CannedCobordismImplData *od = (CannedCobordismImplData *)other->impl_data;
  return CannedCobordismImpl_compose_vertical(sd, od);
}

static CannedCobordism *impl_compose_partial(const CannedCobordism *self,
                                             int start,
                                             const CannedCobordism *cc,
                                             int cstart, int nc) {
  CannedCobordismImplData *sd = (CannedCobordismImplData *)self->impl_data;
  CannedCobordismImplData *cd = (CannedCobordismImplData *)cc->impl_data;
  return CannedCobordismImpl_compose_horizontal(sd, start, cd, cstart, nc);
}

static bool impl_isIsomorphism(const CannedCobordism *self) {
  return CannedCobordismImpl_isIsomorphism(
      (CannedCobordismImplData *)self->impl_data);
}

static void impl_reverseMaps(CannedCobordism *self) {
  CannedCobordismImpl_reverseMaps((CannedCobordismImplData *)self->impl_data);
}

static void impl_free(CannedCobordism *self) {
  CannedCobordismImpl_free((CannedCobordismImplData *)self->impl_data);
  free(self);
}

/* ================================================================
 *  CannedCobordismImpl_as_CannedCobordism
 * ================================================================ */
CannedCobordism *
CannedCobordismImpl_as_CannedCobordism(CannedCobordismImplData *impl) {
  CannedCobordism *cc = (CannedCobordism *)malloc(sizeof(CannedCobordism));
  if (!cc)
    return NULL;
  cc->source = impl->top;
  cc->target = impl->bottom;
  cc->impl_data = impl;
  cc->compose = impl_compose;
  cc->compose_partial = impl_compose_partial;
  cc->isIsomorphism = impl_isIsomorphism;
  cc->reverseMaps = impl_reverseMaps;
  cc->free_impl = impl_free;
  return cc;
}

/* ================================================================
 *  CannedCobordismImpl_create
 *
 *  Purpose:
 *      Constructor — traces boundary components through top/bottom
 *      pairings.
 *
 *  Arguments:
 *      top    — source Cap
 *      bottom — target Cap (must have same n)
 *
 *  Output:
 *      Heap-allocated impl data with component[] filled,
 *      connectedComponent[] allocated but needing external setup.
 *
 *  Method:
 *      For each unvisited edge, trace the cycle
 *      edge -> top.pairings -> bottom.pairings -> ...
 *      marking all visited edges with the boundary component index.
 * ================================================================ */
CannedCobordismImplData *CannedCobordismImpl_create(Cap *top, Cap *bottom) {
  assert(top->n == bottom->n);
  int n = top->n;

  int nbc_max = n + top->ncycles + bottom->ncycles + 2;
  size_t footprint = sizeof(CannedCobordismImplData) +
                     (size_t)n * sizeof(int8_t) +
                     (size_t)nbc_max * sizeof(int8_t);

  char *mem_block = (char *)calloc(1, footprint);
  CannedCobordismImplData *d = (CannedCobordismImplData *)mem_block;
  mem_block += sizeof(CannedCobordismImplData);

  d->component = (int8_t *)mem_block;
  mem_block += (size_t)n * sizeof(int8_t);

  d->connectedComponent = (int8_t *)mem_block;

  d->n = n;
  d->hpower = 0;
  d->top = top;
  d->bottom = bottom;

  memset(d->component, -1, (size_t)n * sizeof(int8_t));

  d->nbc = 0;
  for (int i = 0; i < n; i++) {
    if (d->component[i] == -1) {
      int j = i;
      do {
        d->component[j] = (int8_t)d->nbc;
        j = top->pairings[j];
        d->component[j] = (int8_t)d->nbc;
        j = bottom->pairings[j];
      } while (j != i);
      d->nbc++;
    }
  }

  d->offtop = d->nbc;
  d->nbc += (int8_t)top->ncycles;
  d->offbot = d->nbc;
  d->nbc += (int8_t)bottom->ncycles;

  d->dots = NULL;
  d->genus = NULL;
  d->boundaryComponents = NULL;
  d->bc_sizes = NULL;
  d->edges = NULL;
  d->edge_sizes = NULL;

  d->reverse_maps_done = false;
  d->hashcode = 0;

  return d;
}

/* ================================================================
 *  CannedCobordismImpl_free
 *
 *  Purpose:
 *      Frees a CannedCobordismImplData object and all associated memory.
 *
 *  Method:
 *      1. Frees the lazily-allocated dots and genus arrays.
 *      2. Recursively frees the lazy reverse mapping structures
 * (boundaryComponents, edges).
 *      3. Frees the primary Flat Allocation block (the impl pointer itself),
 *         which implicitly frees the embedded component and connectedComponent
 * arrays.
 * ================================================================ */
void CannedCobordismImpl_free(CannedCobordismImplData *impl) {
  if (!impl)
    return;
  free(impl->dots);
  free(impl->genus);
  if (impl->boundaryComponents) {
    for (int i = 0; i < impl->ncc; i++)
      free(impl->boundaryComponents[i]);
    free(impl->boundaryComponents);
  }
  free(impl->bc_sizes);
  if (impl->edges) {
    for (int i = 0; i < impl->offtop; i++)
      free(impl->edges[i]);
    free(impl->edges);
  }
  free(impl->edge_sizes);
  free(impl);
}

/* ================================================================
 *  CannedCobordismImpl_reverseMaps
 *
 *  Purpose:
 *      Lazily fills boundaryComponents[][] and edges[][].
 *
 *  Arguments:
 *      impl — impl data to fill
 *
 *  Output:
 *      None (modifies impl in place).
 *
 *  Method:
 *      Count BCs per connected component, allocate sub-arrays, fill.
 *      Count edges per mixed BC, allocate sub-arrays, fill.
 * ================================================================ */
void CannedCobordismImpl_reverseMaps(CannedCobordismImplData *impl) {
  if (impl->reverse_maps_done)
    return;

  /* boundaryComponents: which BCs belong to each CC */
  int numBC[impl->ncc];
  memset(numBC, 0, (size_t)impl->ncc * sizeof(int));
  for (int i = 0; i < impl->nbc; i++)
    numBC[impl->connectedComponent[i]]++;

  impl->boundaryComponents =
      (int8_t **)malloc((size_t)impl->ncc * sizeof(int8_t *));
  impl->bc_sizes = (int *)malloc((size_t)impl->ncc * sizeof(int));
  for (int i = 0; i < impl->ncc; i++) {
    impl->bc_sizes[i] = numBC[i];
    impl->boundaryComponents[i] =
        (int8_t *)malloc((size_t)numBC[i] * sizeof(int8_t));
  }

  int j[impl->ncc];
  memset(j, 0, (size_t)impl->ncc * sizeof(int));
  for (int8_t i = 0; i < impl->nbc; i++) {
    int k = impl->connectedComponent[i];
    impl->boundaryComponents[k][j[k]++] = i;
  }

  /* edges: which edges belong to each mixed BC */
  int nedges[impl->n + 1];
  memset(nedges, 0, (size_t)(impl->n + 1) * sizeof(int));
  for (int i = 0; i < impl->n; i++)
    nedges[impl->component[i]]++;

  impl->edges = (int8_t **)malloc((size_t)impl->offtop * sizeof(int8_t *));
  impl->edge_sizes = (int *)malloc((size_t)impl->offtop * sizeof(int));
  for (int i = 0; i < impl->offtop; i++) {
    impl->edge_sizes[i] = nedges[i];
    impl->edges[i] = (int8_t *)malloc((size_t)nedges[i] * sizeof(int8_t));
  }

  int j_edges[impl->offtop];
  memset(j_edges, 0, (size_t)impl->offtop * sizeof(int));
  for (int i = 0; i < impl->n; i++) {
    int k = impl->component[i];
    impl->edges[k][j_edges[k]++] = (int8_t)i;
  }

  impl->reverse_maps_done = true;
}

/* ================================================================
 *  CannedCobordismImpl_check
 * ================================================================ */
bool CannedCobordismImpl_check(const CannedCobordismImplData *impl) {
  if (impl->nbc > 0 && impl->ncc < 1)
    return false;
  if (impl->ncc < 0)
    return false;
  for (int i = 0; i < impl->nbc; i++)
    if (impl->connectedComponent[i] < 0 ||
        impl->connectedComponent[i] >= impl->ncc)
      return false;
  return true;
}

/* ================================================================
 *  CannedCobordismImpl_equals
 * ================================================================ */
bool CannedCobordismImpl_equals(const CannedCobordismImplData *a,
                                const CannedCobordismImplData *b) {
  if (a == b)
    return true;
  if (a->ncc != b->ncc)
    return false;
  if (memcmp(a->dots, b->dots, (size_t)a->ncc) != 0)
    return false;
  if (memcmp(a->genus, b->genus, (size_t)a->ncc) != 0)
    return false;
  if (a->nbc != b->nbc)
    return false;
  if (memcmp(a->connectedComponent, b->connectedComponent, (size_t)a->nbc) != 0)
    return false;
  if (a->n != b->n)
    return false;
  if (a->hpower != b->hpower)
    return false;
  if (memcmp(a->component, b->component, (size_t)a->n) != 0)
    return false;
  if (!Cap_equals(a->top, b->top))
    return false;
  if (!Cap_equals(a->bottom, b->bottom))
    return false;
  return true;
}

/* ================================================================
 *  CannedCobordismImpl_hashCode
 * ================================================================ */
static int arrays_hashCode_i8(const int8_t *arr, int len) {
  int r = 1;
  for (int i = 0; i < len; i++)
    r = 31 * r + arr[i];
  return r;
}

int CannedCobordismImpl_hashCode(CannedCobordismImplData *impl) {
  if (impl->hashcode != 0)
    return impl->hashcode;
  int r = impl->n;
  r += Cap_hashCode(impl->top);
  r += Cap_hashCode(impl->bottom) << 1;
  r += arrays_hashCode_i8(impl->connectedComponent, impl->nbc) << 3;
  r += arrays_hashCode_i8(impl->dots, impl->ncc) << 5;
  r += arrays_hashCode_i8(impl->genus, impl->ncc) << 7;
  r += impl->hpower << 9;
  impl->hashcode = r;
  return r;
}

/* ================================================================
 *  CannedCobordismImpl_compareTo
 * ================================================================ */
int CannedCobordismImpl_compareTo(const CannedCobordismImplData *a,
                                  const CannedCobordismImplData *b) {
  if (a == b)
    return 0;
  if (a->nbc != b->nbc)
    return a->nbc - b->nbc;
  for (int i = 0; i < a->nbc; i++)
    if (a->connectedComponent[i] != b->connectedComponent[i])
      return a->connectedComponent[i] - b->connectedComponent[i];
  if (a->ncc != b->ncc)
    return a->ncc - b->ncc;
  for (int i = 0; i < a->ncc; i++) {
    if (a->dots[i] != b->dots[i])
      return a->dots[i] - b->dots[i];
    if (a->genus[i] != b->genus[i])
      return a->genus[i] - b->genus[i];
  }
  if (a->n != b->n)
    return a->n - b->n;
  if (a->hpower != b->hpower)
    return a->hpower - b->hpower;
  for (int i = 0; i < a->n; i++) {
    if (a->top->pairings[i] != b->top->pairings[i])
      return a->top->pairings[i] - b->top->pairings[i];
    if (a->bottom->pairings[i] != b->bottom->pairings[i])
      return a->bottom->pairings[i] - b->bottom->pairings[i];
  }
  if (a->top->ncycles != b->top->ncycles)
    return a->top->ncycles - b->top->ncycles;
  if (a->bottom->ncycles != b->bottom->ncycles)
    return a->bottom->ncycles - b->bottom->ncycles;
  return 0;
}

/* ================================================================
 *  CannedCobordismImpl_isIsomorphism
 * ================================================================ */
bool CannedCobordismImpl_isIsomorphism(const CannedCobordismImplData *impl) {
  if (!Cap_equals(impl->top, impl->bottom))
    return false;
  if (impl->top->ncycles != 0)
    return false; /* not delooped */
  if (impl->nbc != impl->ncc)
    return false;
  if (impl->hpower != 0)
    return false;
  for (int i = 0; i < impl->nbc; i++) {
    if (impl->connectedComponent[i] != i)
      return false;
    if (impl->dots[i] != 0)
      return false;
    if (impl->genus[i] != 0)
      return false;
  }
  return true;
}

/* ================================================================
 *  CannedCobordismImpl_isomorphism (static factory)
 * ================================================================ */
CannedCobordism *CannedCobordismImpl_isomorphism(Cap *c) {
  assert(c->ncycles == 0);
  CannedCobordismImplData *d = CannedCobordismImpl_create(c, c);
  d->ncc = d->nbc;
  /* connectedComponent[i] = i */
  for (int8_t i = 0; i < d->nbc; i++)
    d->connectedComponent[i] = i;

  d->dots = (int8_t *)calloc((size_t)d->ncc, sizeof(int8_t));
  d->genus = (int8_t *)calloc((size_t)d->ncc, sizeof(int8_t));

  return CannedCobordismImpl_as_CannedCobordism(d);
}

/* ================================================================
 *  Helper: merge connected component labels
 *
 *  Purpose:
 *      Replaces all occurrences of 'old_val' with 'new_val' in
 *      ret->connectedComponent and midConComp arrays.
 *
 *  Arguments:
 *      ret_cc     — ret's connectedComponent array
 *      ret_nbc    — size of ret_cc
 *      midConComp — middle cycle CC array
 *      mid_len    — size of midConComp
 *      old_val    — value to replace
 *      new_val    — replacement value
 *
 *  Output:
 *      None (modifies arrays in place).
 *
 *  Method:
 *      Linear scan of both arrays.
 * ================================================================ */
static void merge_labels(int8_t *ret_cc, int ret_nbc, int *midConComp,
                         int mid_len, int old_val, int new_val) {
  for (int i = 0; i < ret_nbc; i++)
    if (ret_cc[i] == old_val)
      ret_cc[i] = (int8_t)new_val;
  for (int i = 0; i < mid_len; i++)
    if (midConComp[i] == old_val)
      midConComp[i] = new_val;
}

/* ================================================================
 *  CannedCobordismImpl_compose_vertical
 *
 *  Purpose:
 *      Vertical composition: stacks cc on top of this.
 *      this.top must equal cc.bottom.
 *
 *  Arguments:
 *      this_d — bottom cobordism impl data
 *      cc_d   — top cobordism impl data
 *
 *  Output:
 *      New CannedCobordism (vtable-wrapped) representing cc ∘ this.
 *
 *  Method:
 *      1. reverseMaps on both inputs.
 *      2. Create ret with top=cc.top, bottom=this.bottom.
 *      3. Initialize each ret BC as its own CC.
 *      4. For each CC in 'this', find matching ret CC via boundary
 *         components / edges / middle cycles, merge labels.
 *      5. Same for CC in 'cc'.
 *      6. Relabel ret CCs in ascending order.
 *      7. Compute genus via Euler characteristic for each ret CC.
 *      8. Collect unconnected components.
 * ================================================================ */
CannedCobordism *
CannedCobordismImpl_compose_vertical(CannedCobordismImplData *this_d,
                                     CannedCobordismImplData *cc_d) {
  assert(Cap_equals(this_d->top, cc_d->bottom));
  CannedCobordismImpl_reverseMaps(this_d);
  CannedCobordismImpl_reverseMaps(cc_d);

  CannedCobordismImplData *ret =
      CannedCobordismImpl_create(cc_d->top, this_d->bottom);
  ret->hpower = this_d->hpower + cc_d->hpower;

  for (int8_t i = 0; i < ret->nbc; i++)
    ret->connectedComponent[i] = i;

  int mid_len = this_d->top->ncycles;
  int midConComp[mid_len + 1];
  for (int i = 0; i < mid_len; i++)
    midConComp[i] = -2 - i;

  int alloc = this_d->ncc + cc_d->ncc + ret->nbc + 4;
  int8_t rdots[alloc];
  memset(rdots, 0, alloc);
  int8_t mdots[alloc];
  memset(mdots, 0, alloc);
  int8_t udots[alloc];
  memset(udots, 0, alloc);
  int8_t ugenus[alloc];
  memset(ugenus, 0, alloc);
  int unconnected = 0;

  /* --- Merge CCs from 'this' --- */
  for (int i = 0; i < this_d->ncc; i++) {
    int reti = -1;
    for (int j = 0; j < this_d->bc_sizes[i]; j++) {
      int k = this_d->boundaryComponents[i][j];
      if (k < this_d->offtop) {
        for (int l = 0; l < this_d->edge_sizes[k]; l++) {
          reti = ret->connectedComponent[ret->component[(
              int)(this_d->edges[k][l])]];
          if (reti >= 0)
            break;
        }
      } else if (k < this_d->offbot)
        reti = midConComp[k - this_d->offtop];
      else
        reti = ret->connectedComponent[k - this_d->offbot + ret->offbot];
      if (reti >= 0)
        break;
    }

    if (reti == -1) {
      ugenus[unconnected] = this_d->genus[i];
      udots[unconnected++] = this_d->dots[i];
      continue;
    } else if (reti >= 0)
      rdots[reti] += this_d->dots[i];
    else
      mdots[-2 - reti] += this_d->dots[i];

    for (int j = 0; j < this_d->bc_sizes[i]; j++) {
      int bc = this_d->boundaryComponents[i][j];
      if (bc < this_d->offtop) {
        for (int k2 = 0; k2 < this_d->edge_sizes[bc]; k2++) {
          int rt = ret->connectedComponent[ret->component[(
              int)(this_d->edges[bc][k2])]];
          if (rt != reti) {
            merge_labels(ret->connectedComponent, ret->nbc, midConComp, mid_len,
                         rt, reti);
            rdots[reti] += rdots[rt];
          }
        }
      } else if (bc < this_d->offbot) {
        int mtest = midConComp[bc - this_d->offtop];
        if (mtest != reti) {
          if (mtest >= 0)
            for (int k2 = 0; k2 < ret->nbc; k2++)
              if (ret->connectedComponent[k2] == mtest)
                ret->connectedComponent[k2] = (int8_t)reti;
          for (int k2 = 0; k2 < mid_len; k2++)
            if (midConComp[k2] == mtest)
              midConComp[k2] = reti;
          if (reti >= 0) {
            if (mtest >= 0)
              rdots[reti] += rdots[mtest];
            else
              rdots[reti] += mdots[-2 - mtest];
          } else {
            if (mtest >= 0)
              mdots[-2 - reti] += rdots[mtest];
            else
              mdots[-2 - reti] += mdots[-2 - mtest];
          }
        }
      } else {
        int rtest = ret->connectedComponent[bc - this_d->offbot + ret->offbot];
        if (rtest != reti) {
          merge_labels(ret->connectedComponent, ret->nbc, midConComp, mid_len,
                       rtest, reti);
          rdots[reti] += rdots[rtest];
        }
      }
    }
  }

  /* --- Merge CCs from 'cc' (same logic, different offsets) --- */
  for (int i = 0; i < cc_d->ncc; i++) {
    int reti = -1;
    for (int j = 0; j < cc_d->bc_sizes[i]; j++) {
      int k = cc_d->boundaryComponents[i][j];
      if (k < cc_d->offtop) {
        for (int l = 0; l < cc_d->edge_sizes[k]; l++) {
          reti =
              ret->connectedComponent[ret->component[(int)(cc_d->edges[k][l])]];
          if (reti >= 0)
            break;
        }
      } else if (k < cc_d->offbot)
        reti = ret->connectedComponent[k - cc_d->offtop + ret->offtop];
      else
        reti = midConComp[k - cc_d->offbot];
      if (reti >= 0)
        break;
    }

    if (reti == -1) {
      ugenus[unconnected] = cc_d->genus[i];
      udots[unconnected++] = cc_d->dots[i];
      continue;
    } else if (reti >= 0)
      rdots[reti] += cc_d->dots[i];
    else
      mdots[-2 - reti] += cc_d->dots[i];

    for (int j = 0; j < cc_d->bc_sizes[i]; j++) {
      int bc = cc_d->boundaryComponents[i][j];
      if (bc < cc_d->offtop) {
        for (int k2 = 0; k2 < cc_d->edge_sizes[bc]; k2++) {
          int rt = ret->connectedComponent[ret->component[(
              int)(cc_d->edges[bc][k2])]];
          if (rt != reti) {
            merge_labels(ret->connectedComponent, ret->nbc, midConComp, mid_len,
                         rt, reti);
            rdots[reti] += rdots[rt];
          }
        }
      } else if (bc < cc_d->offbot) {
        int rtest = ret->connectedComponent[bc - cc_d->offtop + ret->offtop];
        if (rtest != reti) {
          merge_labels(ret->connectedComponent, ret->nbc, midConComp, mid_len,
                       rtest, reti);
          rdots[reti] += rdots[rtest];
        }
      } else {
        int mtest = midConComp[bc - cc_d->offbot];
        if (mtest != reti) {
          if (mtest >= 0)
            for (int k2 = 0; k2 < ret->nbc; k2++)
              if (ret->connectedComponent[k2] == mtest)
                ret->connectedComponent[k2] = (int8_t)reti;
          for (int k2 = 0; k2 < mid_len; k2++)
            if (midConComp[k2] == mtest)
              midConComp[k2] = reti;
          if (reti >= 0) {
            if (mtest >= 0)
              rdots[reti] += rdots[mtest];
            else
              rdots[reti] += mdots[-2 - mtest];
          } else {
            if (mtest >= 0)
              mdots[-2 - reti] += rdots[mtest];
            else
              mdots[-2 - reti] += mdots[-2 - mtest];
          }
        }
      }
    }
  }

  /* --- Relabel ret CCs in ascending order --- */
  ret->ncc = 0;
  for (int i = 0; i < ret->nbc; i++) {
    if (ret->connectedComponent[i] > ret->ncc) {
      int8_t j = ret->connectedComponent[i];
      for (int k = i; k < ret->nbc; k++) {
        if (ret->connectedComponent[k] == j)
          ret->connectedComponent[k] = ret->ncc;
        else if (ret->connectedComponent[k] == ret->ncc)
          ret->connectedComponent[k] = j;
      }
      for (int k = 0; k < mid_len; k++) {
        if (midConComp[k] == j)
          midConComp[k] = ret->ncc;
        else if (midConComp[k] == ret->ncc)
          midConComp[k] = j;
      }
      int8_t tmp = rdots[ret->ncc];
      rdots[ret->ncc] = rdots[j];
      rdots[j] = tmp;
      ret->ncc++;
    } else if (ret->connectedComponent[i] == ret->ncc)
      ret->ncc++;
  }

  /* --- Genus calculation via Euler characteristic --- */
  CannedCobordismImpl_reverseMaps(ret);
  int8_t rgenus[ret->ncc + 1];
  memset(rgenus, 0, ret->ncc + 1);

  for (int i = 0; i < ret->ncc; i++) {
    int b = ret->bc_sizes[i];
    int x = 0;
    /* Euler contributions from 'this' */
    for (int j = 0; j < this_d->ncc; j++) {
      for (int k = 0; k < this_d->bc_sizes[j]; k++) {
        bool found = false;
        int bc = this_d->boundaryComponents[j][k];
        if (bc < this_d->offtop) {
          for (int l = 0; l < this_d->edge_sizes[bc]; l++)
            if (ret->connectedComponent[ret->component[(
                    int)(this_d->edges[bc][l])]] == i) {
              found = true;
              break;
            }
        } else if (bc < this_d->offbot) {
          if (midConComp[bc - this_d->offtop] == i)
            found = true;
        } else {
          if (ret->connectedComponent[bc - this_d->offbot + ret->offbot] == i)
            found = true;
        }
        if (found) {
          x += 2 - 2 * this_d->genus[j] - this_d->bc_sizes[j];
          int njoins = 0;
          for (int l = 0; l < this_d->bc_sizes[j]; l++) {
            if (this_d->boundaryComponents[j][l] < this_d->offtop)
              njoins += this_d->edge_sizes[this_d->boundaryComponents[j][l]];
            else
              break;
          }
          assert(njoins % 2 == 0);
          x -= njoins / 2;
          break;
        }
      }
    }
    /* Euler contributions from 'cc' */
    for (int j = 0; j < cc_d->ncc; j++) {
      for (int k = 0; k < cc_d->bc_sizes[j]; k++) {
        bool found = false;
        int bc = cc_d->boundaryComponents[j][k];
        if (bc < cc_d->offtop) {
          for (int l = 0; l < cc_d->edge_sizes[bc]; l++)
            if (ret->connectedComponent[ret->component[(
                    int)(cc_d->edges[bc][l])]] == i) {
              found = true;
              break;
            }
        } else if (bc < cc_d->offbot) {
          if (ret->connectedComponent[bc - cc_d->offtop + ret->offtop] == i)
            found = true;
        } else {
          if (midConComp[bc - cc_d->offbot] == i)
            found = true;
        }
        if (found) {
          x += 2 - 2 * cc_d->genus[j] - cc_d->bc_sizes[j];
          break;
        }
      }
    }
    int g = 2 - b - x;
    assert(g % 2 == 0 && g >= 0);
    rgenus[i] = (int8_t)(g / 2);
  }

  /* Handle middle cycles that became unconnected */
  for (int i = 0; i < mid_len; i++) {
    bool found = false;
    int g = 1;
    for (int j = 0; j < mid_len; j++)
      if (midConComp[j] == -2 - i) {
        found = true;
        g++;
      }
    if (found) {
      for (int j = 0; j < this_d->ncc; j++)
        for (int k = this_d->offtop; k < this_d->offbot; k++)
          if (midConComp[k - this_d->offtop] == -2 - i)
            if (this_d->connectedComponent[k] == j) {
              g += this_d->genus[j] - 1;
              break;
            }
      for (int j = 0; j < cc_d->ncc; j++)
        for (int k = cc_d->offbot; k < cc_d->nbc; k++)
          if (midConComp[k - cc_d->offbot] == -2 - i)
            if (cc_d->connectedComponent[k] == j) {
              g += cc_d->genus[j] - 1;
              break;
            }
      ugenus[unconnected] = (int8_t)g;
      udots[unconnected++] = mdots[i];
    }
  }

  /* --- Finalize dots and genus arrays --- */
  int rncc = ret->ncc;
  ret->ncc += (int8_t)unconnected;
  ret->dots = (int8_t *)malloc((size_t)ret->ncc * sizeof(int8_t));
  memcpy(ret->dots, rdots, (size_t)rncc);
  memcpy(ret->dots + rncc, udots, (size_t)unconnected);
  ret->genus = (int8_t *)malloc((size_t)ret->ncc * sizeof(int8_t));
  memcpy(ret->genus, rgenus, (size_t)rncc);
  memcpy(ret->genus + rncc, ugenus, (size_t)unconnected);

  /* Undo reverseMaps if unconnected changed ncc */
  if (unconnected > 0) {
    if (ret->boundaryComponents) {
      for (int i = 0; i < rncc; i++)
        free(ret->boundaryComponents[i]);
      free(ret->boundaryComponents);
      ret->boundaryComponents = NULL;
    }
    free(ret->bc_sizes);
    ret->bc_sizes = NULL;
    if (ret->edges) {
      for (int i = 0; i < ret->offtop; i++)
        free(ret->edges[i]);
      free(ret->edges);
      ret->edges = NULL;
    }
    free(ret->edge_sizes);
    ret->edge_sizes = NULL;
    ret->reverse_maps_done = false;
  }

  assert(CannedCobordismImpl_check(ret));

  return CannedCobordismImpl_as_CannedCobordism(ret);
}

/* ================================================================
 *  CannedCobordismImpl_compose_horizontal
 *
 *  Purpose:
 *      Horizontal composition: composes cc next to this.
 *
 *  Arguments:
 *      this_d — left/this cobordism impl data
 *      start  — offset in this cobordism boundary
 *      cc_d   — right/that cobordism impl data
 *      cstart — offset in cc cobordism boundary
 *      nc     — number of edges to join
 *
 *  Output:
 *      New CannedCobordism (vtable-wrapped) representing cc ∘ this
 * horizontally.
 *
 *  Method:
 *      Uses Cap_compose on top and bottom caps.
 *      Traces boundary components across horizontal connections,
 *      builds a new graph of connected components, Euler characteristic
 * calculations to deduce new genus, and collects unconnected dots/genus.
 * ================================================================ */
CannedCobordism *
CannedCobordismImpl_compose_horizontal(CannedCobordismImplData *this_d,
                                       int start, CannedCobordismImplData *cc_d,
                                       int cstart, int nc) {
  assert(CannedCobordismImpl_check(this_d) && CannedCobordismImpl_check(cc_d));

  int tjoins[nc + 1];
  memset(tjoins, 0, (nc + 1) * sizeof(int));
  int bjoins[nc + 1];
  memset(bjoins, 0, (nc + 1) * sizeof(int));
  Cap *rtop = Cap_compose(this_d->top, start, cc_d->top, cstart, nc, tjoins);
  Cap *rbot =
      Cap_compose(this_d->bottom, start, cc_d->bottom, cstart, nc, bjoins);

  CannedCobordismImplData *ret = CannedCobordismImpl_create(rtop, rbot);
  ret->hpower = this_d->hpower + cc_d->hpower;

  for (int8_t i = 0; i < ret->nbc; i++)
    ret->connectedComponent[i] = i;

  int midConComp[nc + 1];
  for (int i = 0; i < nc; i++)
    midConComp[i] = -2 - i;

  int alloc = this_d->ncc + cc_d->ncc + ret->nbc + nc + 4;
  int8_t rdots[alloc];
  memset(rdots, 0, alloc);
  int8_t mdots[alloc];
  memset(mdots, 0, alloc);
  int8_t udots[alloc];
  memset(udots, 0, alloc);
  int8_t ugenus[alloc];
  memset(ugenus, 0, alloc);
  int unconnected = 0;

  for (int i = 0; i < this_d->ncc; i++) {
    int reti = -1;
    for (int j = this_d->offtop; j < this_d->nbc; j++) {
      if (this_d->connectedComponent[j] == i) {
        int retj;
        if (j < this_d->offbot)
          retj = j - this_d->offtop + ret->offtop;
        else
          retj = j - this_d->offbot + ret->offbot;
        reti = ret->connectedComponent[retj];
        break;
      }
    }
    if (reti == -1) {
      for (int j = 0; j < this_d->n; j++) {
        int thisj = (j + start + nc) % this_d->n;
        if (this_d->connectedComponent[this_d->component[thisj]] == i) {
          if (j < this_d->n - nc)
            reti = ret->connectedComponent[ret->component[j]];
          else
            reti = midConComp[j - this_d->n + nc];
          if (reti >= 0)
            break;
        }
      }
    }

    if (reti >= 0)
      rdots[reti] += this_d->dots[i];
    else if (reti < -1)
      mdots[-2 - reti] += this_d->dots[i];

    if (reti != -1) {
      for (int j = 0; j < this_d->offtop; j++) {
        if (this_d->connectedComponent[j] == i) {
          for (int k = 0; k < this_d->n; k++) {
            if (this_d->component[k] == j) {
              int retk = (k - start - nc + 2 * this_d->n) % this_d->n;
              int tmp;
              if (retk < this_d->n - nc)
                tmp = ret->connectedComponent[ret->component[retk]];
              else
                tmp = midConComp[retk - this_d->n + nc];

              if (tmp != reti) {
                merge_labels(ret->connectedComponent, ret->nbc, midConComp, nc,
                             tmp, reti);
                if (tmp >= 0) {
                  if (reti >= 0)
                    rdots[reti] += rdots[tmp];
                  else
                    mdots[-2 - reti] += rdots[tmp];
                } else {
                  if (reti >= 0)
                    rdots[reti] += mdots[-2 - tmp];
                  else
                    mdots[-2 - reti] += mdots[-2 - tmp];
                }
              }
            }
          }
        }
      }
    }

    if (reti >= 0) {
      for (int j = this_d->offtop; j < this_d->nbc; j++) {
        if (this_d->connectedComponent[j] == i) {
          int retj;
          if (j < this_d->offbot)
            retj = j - this_d->offtop + ret->offtop;
          else
            retj = j - this_d->offbot + ret->offbot;
          if (ret->connectedComponent[retj] != reti) {
            int tmp = ret->connectedComponent[retj];
            merge_labels(ret->connectedComponent, ret->nbc, midConComp, nc, tmp,
                         reti);
            rdots[reti] += rdots[tmp];
          }
        }
      }
    }
    if (reti == -1) {
      ugenus[unconnected] = this_d->genus[i];
      udots[unconnected++] = this_d->dots[i];
    }
  }

  for (int i = 0; i < cc_d->ncc; i++) {
    int reti = -1;
    for (int j = cc_d->offtop; j < cc_d->nbc; j++) {
      if (cc_d->connectedComponent[j] == i) {
        int retj;
        if (j < cc_d->offbot)
          retj =
              j - cc_d->offtop + this_d->offbot - this_d->offtop + ret->offtop;
        else
          retj = j - cc_d->offbot + this_d->nbc - this_d->offbot + ret->offbot;
        reti = ret->connectedComponent[retj];
        break;
      }
    }
    if (reti == -1) {
      for (int j = 0; j < cc_d->n; j++) {
        int thisj = (j + cstart + nc) % cc_d->n;
        if (cc_d->connectedComponent[cc_d->component[thisj]] == i) {
          if (j < cc_d->n - nc)
            reti = ret->connectedComponent[ret->component[j + this_d->n - nc]];
          else
            reti = midConComp[cc_d->n - j - 1];
          if (reti >= 0)
            break;
        }
      }
    }

    if (reti >= 0)
      rdots[reti] += cc_d->dots[i];
    else if (reti < -1)
      mdots[-2 - reti] += cc_d->dots[i];

    if (reti != -1) {
      for (int j = 0; j < cc_d->offtop; j++) {
        if (cc_d->connectedComponent[j] == i) {
          for (int k = 0; k < cc_d->n; k++) {
            if (cc_d->component[k] == j) {
              int retk = (k - cstart - nc + 2 * cc_d->n) % cc_d->n;
              int tmp;
              if (retk < cc_d->n - nc)
                tmp = ret->connectedComponent[ret->component[retk + this_d->n -
                                                             nc]];
              else
                tmp = midConComp[cc_d->n - retk - 1];

              if (tmp != reti) {
                merge_labels(ret->connectedComponent, ret->nbc, midConComp, nc,
                             tmp, reti);
                if (tmp >= 0) {
                  if (reti >= 0)
                    rdots[reti] += rdots[tmp];
                  else
                    mdots[-2 - reti] += rdots[tmp];
                } else {
                  if (reti >= 0)
                    rdots[reti] += mdots[-2 - tmp];
                  else
                    mdots[-2 - reti] += mdots[-2 - tmp];
                }
              }
            }
          }
        }
      }
    }

    if (reti >= 0) {
      for (int j = cc_d->offtop; j < cc_d->nbc; j++) {
        if (cc_d->connectedComponent[j] == i) {
          int retj;
          if (j < cc_d->offbot)
            retj = j - cc_d->offtop + this_d->offbot - this_d->offtop +
                   ret->offtop;
          else
            retj =
                j - cc_d->offbot + this_d->nbc - this_d->offbot + ret->offbot;
          if (ret->connectedComponent[retj] != reti) {
            int tmp = ret->connectedComponent[retj];
            merge_labels(ret->connectedComponent, ret->nbc, midConComp, nc, tmp,
                         reti);
            rdots[reti] += rdots[tmp];
          }
        }
      }
    }
    if (reti == -1) {
      ugenus[unconnected] = cc_d->genus[i];
      udots[unconnected++] = cc_d->dots[i];
    }
  }

  for (int i = 0; i < nc; i++) {
    bool found = false;
    for (int j = 0; j < nc; j++) {
      if (tjoins[j] == i) {
        int reti = ret->offbot - 1 - i;
        if (midConComp[j] >= 0)
          ret->connectedComponent[reti] = midConComp[j];
        else {
          ret->connectedComponent[reti] = (int8_t)(ret->nbc - midConComp[j]);
          rdots[ret->nbc - midConComp[j]] = mdots[-2 - midConComp[j]];
          midConComp[j] = (int)(ret->nbc - midConComp[j]);
        }
        found = true;
        break;
      }
    }
    if (!found)
      break;
  }

  for (int i = 0; i < nc; i++) {
    bool found = false;
    for (int j = 0; j < nc; j++) {
      if (bjoins[j] == i) {
        int reti = ret->nbc - 1 - i;
        if (midConComp[j] >= 0)
          ret->connectedComponent[reti] = midConComp[j];
        else {
          ret->connectedComponent[reti] = (int8_t)(ret->nbc - midConComp[j]);
          rdots[ret->nbc - midConComp[j]] = mdots[-2 - midConComp[j]];
          midConComp[j] = (int)(ret->nbc - midConComp[j]);
        }
        found = true;
        break;
      }
    }
    if (!found)
      break;
  }

  int8_t rgenus[ret->nbc + nc + 4];
  memset(rgenus, 0, ret->nbc + nc + 4);
  for (int i = 0; i < ret->nbc + nc + 2; i++) {
    int b = 0;
    for (int j = 0; j < ret->nbc; j++)
      if (ret->connectedComponent[j] == i)
        b++;
    if (b == 0)
      continue;

    int x = 0;
    for (int j = 0; j < this_d->ncc; j++) {
      for (int k = 0; k < this_d->nbc; k++) {
        if (this_d->connectedComponent[k] == j) {
          bool found = false;
          if (k < this_d->offtop) {
            for (int l = 0; l < this_d->n; l++) {
              if (this_d->component[l] == k) {
                int retl = (l - start - nc + 2 * this_d->n) % this_d->n;
                if (retl < this_d->n - nc) {
                  if (ret->connectedComponent[ret->component[retl]] == i) {
                    found = true;
                    break;
                  }
                } else if (midConComp[retl - this_d->n + nc] == i) {
                  found = true;
                  break;
                }
              }
            }
          } else if (k < this_d->offbot) {
            if (ret->connectedComponent[k - this_d->offtop + ret->offtop] == i)
              found = true;
          } else {
            if (ret->connectedComponent[k - this_d->offbot + ret->offbot] == i)
              found = true;
          }
          if (found) {
            x += 2 - 2 * this_d->genus[j];
            for (int l = 0; l < this_d->nbc; l++)
              if (this_d->connectedComponent[l] == j)
                x--;

            int njoins = 0;
            for (int l = 0; l < nc; l++) {
              if (this_d->connectedComponent[this_d->component[(l + start) %
                                                               this_d->n]] == j)
                njoins++;
            }
            x -= njoins;
            break;
          }
        }
      }
    }

    for (int j = 0; j < cc_d->ncc; j++) {
      for (int k = 0; k < cc_d->nbc; k++) {
        if (cc_d->connectedComponent[k] == j) {
          bool found = false;
          if (k < cc_d->offtop) {
            for (int l = 0; l < cc_d->n; l++) {
              if (cc_d->component[l] == k) {
                int retl = (l - cstart - nc + 2 * cc_d->n) % cc_d->n;
                if (retl < cc_d->n - nc) {
                  if (ret->connectedComponent[ret->component[retl + this_d->n -
                                                             nc]] == i) {
                    found = true;
                    break;
                  }
                } else if (midConComp[cc_d->n - retl - 1] == i) {
                  found = true;
                  break;
                }
              }
            }
          } else if (k < cc_d->offbot) {
            if (ret->connectedComponent[k - cc_d->offtop + this_d->offbot -
                                        this_d->offtop + ret->offtop] == i)
              found = true;
          } else {
            if (ret->connectedComponent[k - cc_d->offbot + this_d->nbc -
                                        this_d->offbot + ret->offbot] == i)
              found = true;
          }
          if (found) {
            x += 2 - 2 * cc_d->genus[j];
            for (int l = 0; l < cc_d->nbc; l++)
              if (cc_d->connectedComponent[l] == j)
                x--;
            break;
          }
        }
      }
    }

    int g = 2 - b - x;
    assert(g % 2 == 0 && g >= 0);
    rgenus[i] = (int8_t)(g / 2);
  }

  ret->ncc = 0;
  for (int i = 0; i < ret->nbc; i++) {
    if (ret->connectedComponent[i] > ret->ncc) {
      int8_t j = ret->connectedComponent[i];
      for (int k = i; k < ret->nbc; k++) {
        if (ret->connectedComponent[k] == j)
          ret->connectedComponent[k] = ret->ncc;
        else if (ret->connectedComponent[k] == ret->ncc)
          ret->connectedComponent[k] = j;
      }
      int8_t tmp = rdots[ret->ncc];
      rdots[ret->ncc] = rdots[j];
      rdots[j] = tmp;
      tmp = rgenus[ret->ncc];
      rgenus[ret->ncc] = rgenus[j];
      rgenus[j] = tmp;
      ret->ncc++;
    } else if (ret->connectedComponent[i] == ret->ncc) {
      ret->ncc++;
    }
  }

  int rncc = ret->ncc;
  ret->ncc += (int8_t)unconnected;
  ret->dots = (int8_t *)malloc((size_t)ret->ncc * sizeof(int8_t));
  memcpy(ret->dots, rdots, (size_t)rncc);
  memcpy(ret->dots + rncc, udots, (size_t)unconnected);

  ret->genus = (int8_t *)malloc((size_t)ret->ncc * sizeof(int8_t));
  memcpy(ret->genus, rgenus, (size_t)rncc);
  memcpy(ret->genus + rncc, ugenus, (size_t)unconnected);

  assert(CannedCobordismImpl_check(ret));

  return CannedCobordismImpl_as_CannedCobordism(ret);
}
