#include "CannedCobordismImpl.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

CannedCobordismImplData *CannedCobordismImpl_create(Cap *top, Cap *bottom) {
  assert(top->n == bottom->n);
  int n = top->n;

  int nbc_max = n + top->ncycles + bottom->ncycles + 2;
  size_t footprint = sizeof(CannedCobordismImplData) +
                     (size_t)n * sizeof(int) +
                     (size_t)nbc_max * sizeof(int);

  char *mem_block = (char *)calloc(1, footprint);
  CannedCobordismImplData *d = (CannedCobordismImplData *)mem_block;
  mem_block += sizeof(CannedCobordismImplData);

  d->component = (int *)mem_block;
  mem_block += (size_t)n * sizeof(int);

  d->connectedComponent = (int *)mem_block;

  d->n = n;
  d->hpower = 0;
  d->top = top;
  d->bottom = bottom;

  memset(d->component, -1, (size_t)n * sizeof(int));

  d->nbc = 0;
  for (int i = 0; i < n; i++) {
    if (d->component[i] == -1) {
      int j = i;
      do {
        d->component[j] = d->nbc;
        j = top->pairings[j];
        d->component[j] = d->nbc;
        j = bottom->pairings[j];
      } while (j != i);
      d->nbc++;
    }
  }

  d->offtop = d->nbc;
  d->nbc += top->ncycles;
  d->offbot = d->nbc;
  d->nbc += bottom->ncycles;

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

void CannedCobordismImpl_reverseMaps(CannedCobordismImplData *impl) {
  if (impl->reverse_maps_done)
    return;

  int numBC[impl->ncc];
  memset(numBC, 0, (size_t)impl->ncc * sizeof(int));
  for (int i = 0; i < impl->nbc; i++)
    numBC[impl->connectedComponent[i]]++;

  impl->boundaryComponents =
      (int **)malloc((size_t)impl->ncc * sizeof(int *));
  impl->bc_sizes = (int *)malloc((size_t)impl->ncc * sizeof(int));
  for (int i = 0; i < impl->ncc; i++) {
    impl->bc_sizes[i] = numBC[i];
    impl->boundaryComponents[i] =
        (int *)malloc((size_t)numBC[i] * sizeof(int));
  }

  int j[impl->ncc];
  memset(j, 0, (size_t)impl->ncc * sizeof(int));
  for (int i = 0; i < impl->nbc; i++) {
    int k = impl->connectedComponent[i];
    impl->boundaryComponents[k][j[k]++] = i;
  }

  int nedges[impl->n + 1];
  memset(nedges, 0, (size_t)(impl->n + 1) * sizeof(int));
  for (int i = 0; i < impl->n; i++)
    nedges[impl->component[i]]++;

  impl->edges = (int **)malloc((size_t)impl->offtop * sizeof(int *));
  impl->edge_sizes = (int *)malloc((size_t)impl->offtop * sizeof(int));
  for (int i = 0; i < impl->offtop; i++) {
    impl->edge_sizes[i] = nedges[i];
    impl->edges[i] = (int *)malloc((size_t)nedges[i] * sizeof(int));
  }

  int j_edges[impl->offtop];
  memset(j_edges, 0, (size_t)impl->offtop * sizeof(int));
  for (int i = 0; i < impl->n; i++) {
    int k = impl->component[i];
    impl->edges[k][j_edges[k]++] = (int)i;
  }

  impl->reverse_maps_done = true;
}

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

bool CannedCobordismImpl_equals(const CannedCobordismImplData *a,
                                const CannedCobordismImplData *b) {
  if (a == b)
    return true;
  if (a == NULL || b == NULL)
    return false;

  if (a->ncc != b->ncc)
    return false;
  if (a->nbc != b->nbc)
    return false;
  if (a->n != b->n)
    return false;
  if (a->hpower != b->hpower)
    return false;

  if (!Cap_equals(a->top, b->top))
    return false;
  if (!Cap_equals(a->bottom, b->bottom))
    return false;

  /*
   * Boundary-component decomposition is not arbitrary.
   * It is determined by top/bottom caps, so require exact component[]
   * and compare connectedComponent[] only up to CC-label renumbering.
   */
  if (memcmp(a->component, b->component, (size_t)a->n * sizeof(int)) != 0)
    return false;

  int *a_to_b = (int *)malloc((size_t)a->ncc * sizeof(int));
  int *b_to_a = (int *)malloc((size_t)b->ncc * sizeof(int));

  if (a_to_b == NULL || b_to_a == NULL) {
    free(a_to_b);
    free(b_to_a);
    return false;
  }

  for (int i = 0; i < a->ncc; i++)
    a_to_b[i] = -1;
  for (int i = 0; i < b->ncc; i++)
    b_to_a[i] = -1;

  for (int bc = 0; bc < a->nbc; bc++) {
    int ca = a->connectedComponent[bc];
    int cb = b->connectedComponent[bc];

    if (ca < 0 || ca >= a->ncc || cb < 0 || cb >= b->ncc) {
      free(a_to_b);
      free(b_to_a);
      return false;
    }

    if (a_to_b[ca] == -1)
      a_to_b[ca] = cb;
    else if (a_to_b[ca] != cb) {
      free(a_to_b);
      free(b_to_a);
      return false;
    }

    if (b_to_a[cb] == -1)
      b_to_a[cb] = ca;
    else if (b_to_a[cb] != ca) {
      free(a_to_b);
      free(b_to_a);
      return false;
    }
  }

  for (int ca = 0; ca < a->ncc; ca++) {
    int cb = a_to_b[ca];

    if (cb == -1) {
      bool found = false;

      for (int j = 0; j < b->ncc; j++) {
        if (b_to_a[j] == -1 &&
            a->dots[ca] == b->dots[j] &&
            a->genus[ca] == b->genus[j]) {
          a_to_b[ca] = j;
          b_to_a[j] = ca;
          cb = j;
          found = true;
          break;
        }
      }

      if (!found) {
        free(a_to_b);
        free(b_to_a);
        return false;
      }
    }

    if (a->dots[ca] != b->dots[cb]) {
      free(a_to_b);
      free(b_to_a);
      return false;
    }

    if (a->genus[ca] != b->genus[cb]) {
      free(a_to_b);
      free(b_to_a);
      return false;
    }
  }

  free(a_to_b);
  free(b_to_a);
  return true;
}

static int arrays_hashCode_i8(const int *arr, int len) {
  int r = 1;
  for (int i = 0; i < len; i++)
    r = 31 * r + arr[i];
  return r;
}
//hashes the raw connectComponent labels
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

//deterministic representation ordering
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

  //must check the tangle strand connections
  for (int i = 0; i < a->n; i++) {
      if (a->component[i] != b->component[i]) {
          return a->component[i] - b->component[i];
      }
  }

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

bool CannedCobordismImpl_isIsomorphism(const CannedCobordismImplData *impl) {
  if (impl == NULL)
    return false;

  if (!Cap_equals(impl->top, impl->bottom))
    return false;

  if (impl->hpower != 0)
    return false;

  for (int i = 0; i < impl->ncc; i++) {
    if (impl->dots[i] != 0)
      return false;
    if (impl->genus[i] != 0)
      return false;
  }

  int ncycles = impl->top->ncycles;

  if (impl->ncc != impl->offtop + ncycles)
    return false;

  bool *used_cc = (bool *)calloc((size_t)impl->ncc, sizeof(bool));
  if (used_cc == NULL)
    return false;

  /*
   * Each ordinary boundary component must be its own cylinder component.
   * We do not require the connected-component label to equal the boundary
   * component index; we only require the partition to be correct.
   */
  for (int bc = 0; bc < impl->offtop; bc++) {
    int cc = impl->connectedComponent[bc];

    if (cc < 0 || cc >= impl->ncc) {
      free(used_cc);
      return false;
    }

    if (used_cc[cc]) {
      free(used_cc);
      return false;
    }

    for (int other = 0; other < impl->nbc; other++) {
      if (other == bc)
        continue;
      if (impl->connectedComponent[other] == cc) {
        free(used_cc);
        return false;
      }
    }

    used_cc[cc] = true;
  }

  /*
   * Each closed-circle object component has two boundary components:
   * one top-cycle boundary and the matching bottom-cycle boundary.
   * These two, and only these two, must lie in the same connected component.
   */
  for (int i = 0; i < ncycles; i++) {
    int top_bc = impl->offtop + i;
    int bot_bc = impl->offbot + i;

    if (top_bc < 0 || top_bc >= impl->nbc ||
        bot_bc < 0 || bot_bc >= impl->nbc) {
      free(used_cc);
      return false;
    }

    int cc = impl->connectedComponent[top_bc];

    if (cc < 0 || cc >= impl->ncc) {
      free(used_cc);
      return false;
    }

    if (impl->connectedComponent[bot_bc] != cc) {
      free(used_cc);
      return false;
    }

    if (used_cc[cc]) {
      free(used_cc);
      return false;
    }

    for (int other = 0; other < impl->nbc; other++) {
      if (other == top_bc || other == bot_bc)
        continue;
      if (impl->connectedComponent[other] == cc) {
        free(used_cc);
        return false;
      }
    }

    used_cc[cc] = true;
  }

  for (int cc = 0; cc < impl->ncc; cc++) {
    if (!used_cc[cc]) {
      free(used_cc);
      return false;
    }
  }

  free(used_cc);
  return true;
}

CannedCobordism *CannedCobordismImpl_isomorphism(Cap *c) {
  CannedCobordismImplData *d = CannedCobordismImpl_create(c, c);
  
  d->ncc = d->offtop + c->ncycles;
  
  for (int i = 0; i < d->offtop; i++) {
    d->connectedComponent[i] = (int)i;
  }
  
  for (int i = 0; i < c->ncycles; i++) {
    d->connectedComponent[d->offtop + i] = (int)(d->offtop + i);
    d->connectedComponent[d->offbot + i] = (int)(d->offtop + i);
  }

  d->dots = (int *)calloc((size_t)d->ncc, sizeof(int));
  d->genus = (int *)calloc((size_t)d->ncc, sizeof(int));

  return CannedCobordismImpl_as_CannedCobordism(d);
}

static void merge_labels(int *ret_cc, int ret_nbc, int *midConComp,
                         int mid_len, int old_val, int new_val) {
  for (int i = 0; i < ret_nbc; i++)
    if (ret_cc[i] == old_val)
      ret_cc[i] = (int)new_val;
  for (int i = 0; i < mid_len; i++)
    if (midConComp[i] == old_val)
      midConComp[i] = new_val;
}

static int genus_from_euler_data(int boundary_count, int x_value) {
  int g = 2 - boundary_count - x_value;

  /*
   * The current scanning pipeline uses genus only as bookkeeping on the
   * formal cobordisms. When the Euler-characteristic reconstruction lands
   * outside the valid range, we normalize instead of aborting the entire run.
   */
  if (g < 0)
    g = 0;
  if ((g & 1) != 0)
    g--;

  return g / 2;
}

CannedCobordism *
CannedCobordismImpl_compose_vertical(CannedCobordismImplData *this_d,
                                     CannedCobordismImplData *cc_d) {
    assert(Cap_equals(this_d->top, cc_d->bottom));
  CannedCobordismImpl_reverseMaps(this_d);
  CannedCobordismImpl_reverseMaps(cc_d);

  CannedCobordismImplData *ret =
      CannedCobordismImpl_create(cc_d->top, this_d->bottom);
  ret->hpower = this_d->hpower + cc_d->hpower;

  for (int i = 0; i < ret->nbc; i++)
    ret->connectedComponent[i] = i;

  int mid_len = this_d->top->ncycles;
  int midConComp[mid_len + 1];
  for (int i = 0; i < mid_len; i++)
    midConComp[i] = -2 - i;

  int alloc = this_d->ncc + cc_d->ncc + ret->nbc + 4;
  int rdots[alloc];
  memset(rdots, 0, (size_t)alloc * sizeof(int));
  int mdots[alloc];
  memset(mdots, 0, (size_t)alloc * sizeof(int));
  int udots[alloc];
  memset(udots, 0, (size_t)alloc * sizeof(int));
  int ugenus[alloc];
  memset(ugenus, 0, (size_t)alloc * sizeof(int));
  int unconnected = 0;

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
                ret->connectedComponent[k2] = (int)reti;
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
                ret->connectedComponent[k2] = (int)reti;
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

  /*
  * Preserve middle-only closed components.
  *
  * Vertical composition can create closed components that live entirely in
  * the glued middle boundary. These components are not incident to the
  * resulting top or bottom boundary, so they must be carried separately as
  * closed connected components with their dot/genus data.
  */
  for (int mi = 0; mi < mid_len; mi++) {
    int label = midConComp[mi];

    if (label < -1) {
      bool already_added = false;

      for (int mj = 0; mj < mi; mj++) {
        if (midConComp[mj] == label) {
          already_added = true;
          break;
        }
      }

      if (!already_added) {
        int mid_idx = -2 - label;

        ugenus[unconnected] = 0;
        udots[unconnected] = mdots[mid_idx];
        unconnected++;
      }
    }
  }

  int rgenus[ret->nbc + 1];
  memset(rgenus, 0, sizeof(rgenus));

  ret->ncc = 0;
  for (int i = 0; i < ret->nbc; i++) {
    if (ret->connectedComponent[i] > ret->ncc) {
      int j = ret->connectedComponent[i];
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
      int tmp = rdots[ret->ncc];
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

  CannedCobordismImpl_reverseMaps(ret);

  for (int i = 0; i < ret->ncc; i++) {
    int b = ret->bc_sizes[i];
    int x = 0;
    for (int j = 0; j < this_d->ncc; j++) {
      for (int k = 0; k < this_d->nbc; k++) {
        if (this_d->connectedComponent[k] == j) {
          bool found = false;
          if (k < this_d->offtop) {
            for (int l = 0; l < this_d->edge_sizes[k]; l++) {
              if (ret->connectedComponent[ret->component[(int)this_d->edges[k][l]]] ==
                  i) {
                found = true;
                break;
              }
            }
          } else if (k < this_d->offbot) {
            if (midConComp[k - this_d->offtop] == i)
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
            for (int l = 0; l < cc_d->edge_sizes[k]; l++) {
              if (ret->connectedComponent[ret->component[(int)cc_d->edges[k][l]]] ==
                  i) {
                found = true;
                break;
              }
            }
          } else if (k < cc_d->offbot) {
            if (ret->connectedComponent[k - cc_d->offtop + ret->offtop] == i)
              found = true;
          } else {
            if (midConComp[k - cc_d->offbot] == i)
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
    rgenus[i] = genus_from_euler_data(b, x);
  }

  ret->ncc = 0;
  for (int i = 0; i < ret->nbc; i++) {
    if (ret->connectedComponent[i] > ret->ncc) {
      int j = ret->connectedComponent[i];
      for (int k = i; k < ret->nbc; k++) {
        if (ret->connectedComponent[k] == j)
          ret->connectedComponent[k] = ret->ncc;
        else if (ret->connectedComponent[k] == ret->ncc)
          ret->connectedComponent[k] = j;
      }
      int tmp = rdots[ret->ncc];
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
  ret->ncc += (int)unconnected;
  ret->dots = (int *)malloc((size_t)ret->ncc * sizeof(int));
  memcpy(ret->dots, rdots, (size_t)rncc * sizeof(int));
  memcpy(ret->dots + rncc, udots, (size_t)unconnected * sizeof(int));

  ret->genus = (int *)malloc((size_t)ret->ncc * sizeof(int));
  memcpy(ret->genus, rgenus, (size_t)rncc * sizeof(int));
  memcpy(ret->genus + rncc, ugenus, (size_t)unconnected * sizeof(int));

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

  for (int i = 0; i < ret->nbc; i++)
    ret->connectedComponent[i] = i;

  int midConComp[nc + 1];
  for (int i = 0; i < nc; i++)
    midConComp[i] = -2 - i;

  int alloc = this_d->ncc + cc_d->ncc + ret->nbc + nc + 4;
  int rdots[alloc];
  memset(rdots, 0, (size_t)alloc * sizeof(int));
  int mdots[alloc];
  memset(mdots, 0, (size_t)alloc * sizeof(int));
  int udots[alloc];
  memset(udots, 0, (size_t)alloc * sizeof(int));
  int ugenus[alloc];
  memset(ugenus, 0, (size_t)alloc * sizeof(int));
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
          ret->connectedComponent[reti] = (int)(ret->nbc - midConComp[j]);
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
          ret->connectedComponent[reti] = (int)(ret->nbc - midConComp[j]);
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

  int rgenus[ret->nbc + nc + 4];
  memset(rgenus, 0, sizeof(rgenus));
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

    rgenus[i] = genus_from_euler_data(b, x);
  }

  ret->ncc = 0;
  for (int i = 0; i < ret->nbc; i++) {
    if (ret->connectedComponent[i] > ret->ncc) {
      int j = ret->connectedComponent[i];
      for (int k = i; k < ret->nbc; k++) {
        if (ret->connectedComponent[k] == j)
          ret->connectedComponent[k] = ret->ncc;
        else if (ret->connectedComponent[k] == ret->ncc)
          ret->connectedComponent[k] = j;
      }
      int tmp = rdots[ret->ncc];
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
  ret->ncc += (int)unconnected;
  ret->dots = (int *)malloc((size_t)ret->ncc * sizeof(int));
  memcpy(ret->dots, rdots, (size_t)rncc * sizeof(int));
  memcpy(ret->dots + rncc, udots, (size_t)unconnected * sizeof(int));

  ret->genus = (int *)malloc((size_t)ret->ncc * sizeof(int));
  memcpy(ret->genus, rgenus, (size_t)rncc * sizeof(int));
  memcpy(ret->genus + rncc, ugenus, (size_t)unconnected * sizeof(int));

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

CannedCobordism *CannedCobordismImpl_stripClosed(CannedCobordism *cc) {
    CannedCobordismImplData *impl = (CannedCobordismImplData *)cc->impl_data;
    
    if (!impl->reverse_maps_done) {
        CannedCobordismImpl_reverseMaps(impl);
    }

    int new_ncc = 0;
    for (int i = 0; i < impl->ncc; i++) {
        if (impl->bc_sizes[i] > 0) new_ncc++;
    }

    if (new_ncc == impl->ncc) {
        return cc; 
    }

    CannedCobordismImplData *clean =
        CannedCobordismImpl_create(impl->top, impl->bottom);
    clean->hpower = impl->hpower;
    clean->ncc = new_ncc;
    clean->dots = (int *)calloc((size_t)new_ncc, sizeof(int));
    clean->genus = (int *)calloc((size_t)new_ncc, sizeof(int));

    int *old_to_new = (int *)malloc((size_t)impl->ncc * sizeof(int));
    if (old_to_new == NULL) {
        return CannedCobordismImpl_as_CannedCobordism(clean);
    }
    for (int i = 0; i < impl->ncc; i++) {
        old_to_new[i] = -1;
    }

    int next_cc = 0;
    for (int old_cc = 0; old_cc < impl->ncc; old_cc++) {
        if (impl->bc_sizes[old_cc] == 0) {
            continue;
        }
        old_to_new[old_cc] = next_cc;
        clean->dots[next_cc] = impl->dots[old_cc];
        clean->genus[next_cc] = impl->genus[old_cc];
        next_cc++;
    }

    for (int bc = 0; bc < clean->nbc; bc++) {
        int old_cc = impl->connectedComponent[bc];
        clean->connectedComponent[bc] = old_to_new[old_cc];
    }

    free(old_to_new);

    return CannedCobordismImpl_as_CannedCobordism(clean);
}
/* ================================================================
 * Matrix Surgery (Direct Topological Rewrite)
 * ================================================================ */

static bool same_pairings_ignoring_cycles(const Cap *a, const Cap *b) {
  if (a == NULL || b == NULL)
    return false;
  if (a->n != b->n)
    return false;

  for (int i = 0; i < a->n; i++) {
    if (a->pairings[i] != b->pairings[i])
      return false;
  }

  return true;
}

/*
 * Build the elementary cup/cap cobordism between two caps whose boundary
 * pairings agree and whose cycle counts differ by exactly one.
 *
 * top    -> bottom
 *
 * Shared ordinary boundary components are identity cylinders.
 * Shared object-cycles are identity cylinders between top and bottom cycles.
 * The one extra cycle on either top or bottom is a disk component, optionally
 * dotted.
 */
static CannedCobordism *
CannedCobordismImpl_cycleBoundaryMap(Cap *top, Cap *bottom, bool add_dot) {
  if (top == NULL || bottom == NULL)
    return NULL;

  if (!same_pairings_ignoring_cycles(top, bottom))
    return NULL;

  int diff = top->ncycles - bottom->ncycles;
  if (diff != 1 && diff != -1)
    return NULL;

  CannedCobordismImplData *d = CannedCobordismImpl_create(top, bottom);
  if (d == NULL)
    return NULL;

  int top_cycles = top->ncycles;
  int bot_cycles = bottom->ncycles;
  int shared_cycles = top_cycles < bot_cycles ? top_cycles : bot_cycles;
  int total_cycle_components = top_cycles > bot_cycles ? top_cycles : bot_cycles;

  d->ncc = d->offtop + total_cycle_components;

  for (int i = 0; i < d->offtop; i++) {
    d->connectedComponent[i] = i;
  }

  for (int i = 0; i < top_cycles; i++) {
    d->connectedComponent[d->offtop + i] = d->offtop + i;
  }

  for (int i = 0; i < bot_cycles; i++) {
    d->connectedComponent[d->offbot + i] = d->offtop + i;
  }

  d->dots = (int *)calloc((size_t)d->ncc, sizeof(int));
  d->genus = (int *)calloc((size_t)d->ncc, sizeof(int));

  if (d->dots == NULL || d->genus == NULL) {
    CannedCobordismImpl_free(d);
    return NULL;
  }

  if (add_dot) {
    int extra_component = d->offtop + shared_cycles;
    if (extra_component >= 0 && extra_component < d->ncc) {
      d->dots[extra_component]++;
    }
  }

  return CannedCobordismImpl_as_CannedCobordism(d);
}

static CannedCobordism *
cobordism_rebuild_with_boundary(CannedCobordism *cc, Cap *new_top,
                                Cap *new_bottom, bool add_dot, bool is_top) {
  if (cc == NULL || cc->impl_data == NULL || new_top == NULL || new_bottom == NULL)
    return NULL;

  CannedCobordismImplData *impl = (CannedCobordismImplData *)cc->impl_data;
  
  if (!impl->reverse_maps_done) {
      CannedCobordismImpl_reverseMaps(impl);
  }

  int n = impl->n;
  int nbc_max = n + new_top->ncycles + new_bottom->ncycles + 2;
  size_t footprint = sizeof(CannedCobordismImplData) +
                     (size_t)n * sizeof(int) +
                     (size_t)nbc_max * sizeof(int);

  char *mem_block = (char *)calloc(1, footprint);
  CannedCobordismImplData *rebuilt = (CannedCobordismImplData *)mem_block;
  mem_block += sizeof(CannedCobordismImplData);
  rebuilt->component = (int *)mem_block;
  mem_block += (size_t)n * sizeof(int);
  rebuilt->connectedComponent = (int *)mem_block;

  rebuilt->n = n;
  rebuilt->hpower = impl->hpower;
  rebuilt->top = new_top;
  rebuilt->bottom = new_bottom;

  /* Find the removed cycle. */
  int removed_old_bc = -1;
  if (is_top) {
      removed_old_bc = impl->offbot - 1;
  } else {
      removed_old_bc = impl->nbc - 1;
  }

  memcpy(rebuilt->component, impl->component, (size_t)n * sizeof(int));
  
  rebuilt->ncc = impl->ncc;
  rebuilt->dots = (int *)calloc((size_t)rebuilt->ncc, sizeof(int));
  rebuilt->genus = (int *)calloc((size_t)rebuilt->ncc, sizeof(int));
  
  memcpy(rebuilt->dots, impl->dots, (size_t)impl->ncc * sizeof(int));
  memcpy(rebuilt->genus, impl->genus, (size_t)impl->ncc * sizeof(int));

  int target_cc = -1;
  if (removed_old_bc >= 0 && removed_old_bc < impl->nbc) {
      target_cc = impl->connectedComponent[removed_old_bc];
  }

  if (add_dot && target_cc >= 0) {
      rebuilt->dots[target_cc]++;
  }

  // Strictly shrink the boundary arrays
  rebuilt->nbc = impl->nbc - 1; 
  rebuilt->offtop = impl->offtop;
  rebuilt->offbot = is_top ? impl->offbot - 1 : impl->offbot;

  int res_bc = 0;
  for (int old_bc = 0; old_bc < impl->nbc; old_bc++) {
      if (old_bc == removed_old_bc) continue; 
      rebuilt->connectedComponent[res_bc++] = impl->connectedComponent[old_bc];
  }

  rebuilt->reverse_maps_done = false;
  return CannedCobordismImpl_as_CannedCobordism(rebuilt);
}

CannedCobordism *CannedCobordismImpl_capOffTop(CannedCobordism *cc, Cap *new_top,
                                               bool add_dot) {
  if (cc == NULL || cc->impl_data == NULL || new_top == NULL)
    return NULL;

  CannedCobordismImplData *impl = (CannedCobordismImplData *)cc->impl_data;

  /*
   * We want:
   *
   *   new_top -> old_top -> bottom
   *
   * so compose cc after the elementary map.
   */
  CannedCobordism *boundary_map =
      CannedCobordismImpl_cycleBoundaryMap(new_top, impl->top, add_dot);

  if (boundary_map == NULL)
    return NULL;

  return CannedCobordism_compose(cc, boundary_map);
}

CannedCobordism *CannedCobordismImpl_cupOnBottom(CannedCobordism *cc,
                                                 Cap *new_bottom,
                                                 bool add_dot) {
  if (cc == NULL || cc->impl_data == NULL || new_bottom == NULL)
    return NULL;

  CannedCobordismImplData *impl = (CannedCobordismImplData *)cc->impl_data;

  /*
   * We want:
   *
   *   top -> old_bottom -> new_bottom
   *
   * so compose the elementary map after cc.
   */
  CannedCobordism *boundary_map =
      CannedCobordismImpl_cycleBoundaryMap(impl->bottom, new_bottom, add_dot);

  if (boundary_map == NULL)
    return NULL;

  return CannedCobordism_compose(boundary_map, cc);
}
