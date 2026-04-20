/*
 * Purpose:
 *      Represents a "cap" (planar matching / smoothing) in the
 *      Khovanov homology computation. A cap consists of n boundary
 *      points connected by planar pairings, plus some number of
 *      interior cycles.
 */

#ifndef CAP_H
#define CAP_H

#include <stdbool.h>
#include <stdint.h>

typedef struct Cap {
  int n;         /* Number of boundary edges/points                     */
  int ncycles;   /* Number of interior cycles                           */
  int *pairings; /* Array of size n: pairings[i] = partner of edge i    */
} Cap;

/*
 * Purpose:
 *      Allocates and initialises a new Cap with the given number of
 *      boundary edges and interior cycles. The pairings array is
 *      allocated but left uninitialised — the caller must fill it.
 * Arguments:
 *      n, ncycles
 * Argument descriptions:
 *      n       — number of boundary edges
 *      ncycles — number of interior (closed) cycles
 * Output:
 *      Cap pointer
 * Output description:
 *      Returns a heap-allocated Cap pointer. NULL on allocation failure.
 * Method:
 *      malloc the struct and pairings array; set n and ncycles.
 */
Cap *Cap_create(int n, int ncycles);

/*
 * Purpose:
 *      Frees a Cap and its pairings array.
 * Arguments:
 *      cap
 * Argument descriptions:
 *      cap — the Cap to free (may be NULL)
 * Output:
 *      None
 * Output description:
 *      None
 * Method:
 *      Free pairings, then free the struct.
 */
void Cap_free(Cap *cap);

/*
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
 *      Returns true if a and b have the same n, ncycles, and pairings.
 * Method:
 *      Compare n, ncycles, then memcmp the pairings arrays.
 */
bool Cap_equals(const Cap *a, const Cap *b);

/*
 * Purpose:
 *      Computes a hash code for a Cap.
 * Arguments:
 *      cap
 * Argument descriptions:
 *      cap — the Cap to hash
 * Output:
 *      integer
 * Output description:
 *      Returns an integer hash code.
 * Method:
 *      Iterate pairings with result = 31*result + element for each element.
 * Then add ncycles.
 */
int Cap_hashCode(const Cap *cap);

/*
 * Purpose:
 *      Provides a total ordering on Caps.
 * Arguments:
 *      a, b
 * Argument descriptions:
 *      a — first Cap
 *      b — second Cap
 * Output:
 *      integer
 * Output description:
 *      Returns negative if a < b, 0 if equal, positive if a > b.
 * Method:
 *      Compare n, then ncycles, then element-wise pairings.
 */
int Cap_compareTo(const Cap *a, const Cap *b);

/*
 * Purpose:
 *      Horizontal composition of two Caps, joining nc consecutive
 *      edges starting at position 'start' in this Cap with nc edges
 *      starting at 'cstart' in the other Cap.
 * Arguments:
 *      a, start, b, cstart, nc, joins
 * Argument descriptions:
 *      a      — first Cap (this)
 *      start  — starting edge index in a
 *      b      — second Cap
 *      cstart — starting edge index in b
 *      nc     — number of edges to join
 *      joins  — output array of size nc; filled with information
 *               about newly created cycles. joins[i] = -1 if edge i
 *               was consumed by an existing path; otherwise joins[i]
 *               is the index of the new cycle it became part of.
 *               May be NULL if join info is not needed.
 * Output:
 *      Cap pointer
 * Output description:
 *      Returns a newly allocated Cap representing the composition.
 * Method:
 *      Remap edges from both caps into a combined pairing array of size
 *      a.n + b.n - 2*nc, then resolve all cross-connections through the
 *      join edges, and finally detect any newly created cycles.
 */
Cap *Cap_compose(const Cap *a, int start, const Cap *b, int cstart, int nc,
                 int *joins);

#endif /* CAP_H */
