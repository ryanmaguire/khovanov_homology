/*
 * Purpose:
 *      Implements the Cap (planar matching / smoothing) used as the
 *      source and target of cobordisms in the Khovanov chain complex.
 */

#include "Cap.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/* ----------------------------------------------------------------
 * Purpose:
 *      Allocates a Cap with n boundary edges and the given cycle count.
 * Arguments:
 *      n, ncycles
 * Argument descriptions:
 *      n       — number of boundary edges
 *      ncycles — number of interior cycles
 * Output:
 *      Cap pointer
 * Output description:
 *      Heap-allocated Cap with uninitialised pairings. NULL on failure.
 * Method:
 *      malloc struct, malloc pairings array of size n.
 * ---------------------------------------------------------------- */
Cap *Cap_create(int n, int ncycles) {
  Cap *cap = (Cap *)malloc(sizeof(Cap));
  if (cap == NULL)
    return NULL;

  cap->n = n;
  cap->ncycles = ncycles;
  cap->pairings = (int *)malloc((size_t)n * sizeof(int));
  if (cap->pairings == NULL && n > 0) {
    free(cap);
    return NULL;
  }
  return cap;
}

/* ----------------------------------------------------------------
 * Purpose:
 *      Frees a Cap and its pairings array.
 * Arguments:
 *      cap
 * Argument descriptions:
 *      cap — Cap to free (may be NULL)
 * Output:
 *      None
 * Output description:
 *      None
 * Method:
 *      Free pairings, then the struct.
 * ---------------------------------------------------------------- */
void Cap_free(Cap *cap) {
  if (cap == NULL)
    return;
  free(cap->pairings);
  free(cap);
}

/* ----------------------------------------------------------------
 * Purpose:
 *      Tests structural equality of two Caps.
 * Arguments:
 *      a, b
 * Argument descriptions:
 *      a — first Cap
 *      b — second Cap
 * Output:
 *      boolean
 * Output description:
 *      true if a and b have identical n, ncycles, and pairings.
 * Method:
 *      Compare n, ncycles, then memcmp the pairings arrays.
 * ---------------------------------------------------------------- */
bool Cap_equals(const Cap *a, const Cap *b) {
  if (a == b)
    return true;
  if (a == NULL || b == NULL)
    return false;
  if (a->n != b->n)
    return false;
  if (a->ncycles != b->ncycles)
    return false;
  if (memcmp(a->pairings, b->pairings, (size_t)a->n * sizeof(int)) != 0)
    return false;
  return true;
}

/* ----------------------------------------------------------------
 * Purpose:
 *      Computes a hash code matching a standard array hash code
 *      of the pairings plus the number of ncycles.
 * Arguments:
 *      cap
 * Argument descriptions:
 *      cap — the Cap to hash
 * Output:
 *      integer
 * Output description:
 *      Integer hash code.
 * Method:
 *      Iterate pairings with result = 31 * result + pairings[i],
 *      starting from result = 1. Then add ncycles.
 * ---------------------------------------------------------------- */
int Cap_hashCode(const Cap *cap) {
  int result = 1;
  for (int i = 0; i < cap->n; i++)
    result = 31 * result + cap->pairings[i];
  return result + cap->ncycles;
}

/* ----------------------------------------------------------------
 * Purpose:
 *      Total ordering on Caps.
 * Arguments:
 *      a, b
 * Argument descriptions:
 *      a — first Cap
 *      b — second Cap
 * Output:
 *      integer
 * Output description:
 *      Negative if a < b, 0 if equal, positive if a > b.
 * Method:
 *      Compare n, then ncycles, then element-wise pairings.
 * ---------------------------------------------------------------- */
int Cap_compareTo(const Cap *a, const Cap *b) {
  if (a->n != b->n)
    return a->n - b->n;
  if (a->ncycles != b->ncycles)
    return a->ncycles - b->ncycles;
  for (int i = 0; i < a->n; i++)
    if (a->pairings[i] != b->pairings[i])
      return a->pairings[i] - b->pairings[i];
  return 0;
}

/* ----------------------------------------------------------------
 * Purpose:
 *      Horizontal composition of two Caps. Joins nc consecutive
 *      edges from each cap and produces a new combined Cap.
 * Arguments:
 *      a, start, b, cstart, nc, joins
 * Argument descriptions:
 *      a      — first Cap
 *      start  — starting edge index in a
 *      b      — second Cap
 *      cstart — starting edge index in b
 *      nc     — number of edges to join
 *      joins  — output array of size nc (may be NULL). Filled with
 *               information about newly created cycles:
 *               joins[i] = -1 if edge i was consumed by an existing
 *               path; otherwise joins[i] = index of the new cycle.
 * Output:
 *      Cap pointer
 * Output description:
 *      A newly heap-allocated Cap representing the composition.
 * Method:
 *      1. Build new pairings of size a.n + b.n - 2*nc by remapping
 *         edges from both caps, marking cross-references as negative.
 *      2. Compute thisjoins[] and cjoins[] — where the joining edges
 *         connect to.
 *      3. Resolve all negative (pending) pairings by following
 *         chains through thisjoins/cjoins.
 *      4. Detect newly created cycles from unconsumed join edges.
 * ---------------------------------------------------------------- */
Cap *Cap_compose(const Cap *a, int start, const Cap *b, int cstart, int nc,
                 int *joins) {
  int an = a->n;
  int bn = b->n;
  int rn = an + bn - 2 * nc;

  Cap *ret = Cap_create(rn, a->ncycles + b->ncycles);

  /* Step 1: Remap edges from cap a (skipping the nc join edges) */
  for (int i = 0; i < an - nc; i++) {
    int ii = (i + start + nc) % an;
    ret->pairings[i] = (a->pairings[ii] - start - nc + 2 * an) % an;
    /* If it pairs with a joining edge, mark as negative */
    if (ret->pairings[i] >= an - nc)
      ret->pairings[i] = -1 - (ret->pairings[i] - an + nc);
  }

  /* Compute thisjoins: where each of a's join edges connects to */
  int *thisjoins = (int *)malloc((size_t)nc * sizeof(int));
  for (int i = 0; i < nc; i++) {
    int ii = (i + start) % an;
    thisjoins[i] = (a->pairings[ii] - start - nc + 2 * an) % an;
    if (thisjoins[i] >= an - nc)
      thisjoins[i] = -1 - (thisjoins[i] - an + nc);
  }

  /* Step 2: Remap edges from cap b (skipping the nc join edges) */
  for (int i = 0; i < bn - nc; i++) {
    int ii = (i + cstart + nc) % bn;
    int j = i + an - nc;
    ret->pairings[j] = (b->pairings[ii] - cstart - nc + 2 * bn) % bn;
    if (ret->pairings[j] >= bn - nc)
      ret->pairings[j] = -1 - (bn - ret->pairings[j] - 1);
    else
      ret->pairings[j] += an - nc;
  }

  /* Compute cjoins: where each of b's join edges connects to */
  int *cjoins = (int *)malloc((size_t)nc * sizeof(int));
  for (int i = 0; i < nc; i++) {
    int ii = (cstart + nc - 1 - i + bn) % bn;
    cjoins[i] = (b->pairings[ii] - cstart - nc + 2 * bn) % bn;
    if (cjoins[i] >= bn - nc)
      cjoins[i] = -1 - (bn - cjoins[i] - 1);
    else
      cjoins[i] += an - nc;
  }

  /* Step 3: Resolve pending (negative) pairings from cap a's side */
  bool *joinsdone = (bool *)calloc((size_t)nc, sizeof(bool));

  for (int i = 0; i < an - nc; i++) {
    if (ret->pairings[i] < 0) {
      while (ret->pairings[i] < 0) {
        int j = -1 - ret->pairings[i];
        joinsdone[j] = true;
        ret->pairings[i] = cjoins[j];
        if (ret->pairings[i] < 0) {
          j = -1 - ret->pairings[i];
          joinsdone[j] = true;
          ret->pairings[i] = thisjoins[j];
        }
      }
      ret->pairings[ret->pairings[i]] = i;
    }
  }

  /* Resolve pending pairings from cap b's side */
  for (int i = an - nc; i < rn; i++) {
    if (ret->pairings[i] < 0) {
      while (ret->pairings[i] < 0) {
        int j = -1 - ret->pairings[i];
        joinsdone[j] = true;
        ret->pairings[i] = thisjoins[j];
        if (ret->pairings[i] < 0) {
          j = -1 - ret->pairings[i];
          joinsdone[j] = true;
          ret->pairings[i] = cjoins[j];
        }
      }
      ret->pairings[ret->pairings[i]] = i;
    }
  }

  /* Step 4: Mark consumed joins and detect newly created cycles */
  int *local_joins = NULL;
  if (joins != NULL) {
    local_joins = joins;
  } else {
    local_joins = (int *)malloc((size_t)nc * sizeof(int));
  }

  for (int i = 0; i < nc; i++) {
    if (joinsdone[i])
      local_joins[i] = -1;
  }

  int nnew = 0;
  for (int i = 0; i < nc; i++) {
    if (!joinsdone[i]) {
      int j = i;
      do {
        joinsdone[j] = true;
        local_joins[j] = nnew;
        j = -1 - thisjoins[j];
        joinsdone[j] = true;
        local_joins[j] = nnew;
        j = -1 - cjoins[j];
      } while (j != i);
      ret->ncycles++;
      nnew++;
    }
  }

  /* Cleanup */
  if (joins == NULL)
    free(local_joins);
  free(thisjoins);
  free(cjoins);
  free(joinsdone);

  return ret;
}
