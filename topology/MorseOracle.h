/*
 *  MorseOracle.h
 *
 *  Topological Discrete Morse Oracle for Fixed-Strand FPT Khovanov Homology.
 *
 *  Purpose:
 *      Implements Phase 1 of the FPT algorithm: given a braid word,
 *      constructs the augmented rhomboid graph G(b), computes the
 *      discrete Morse collapse of its independence complex I(G(b)),
 *      and produces a deterministic schedule of (chain_idx, source_col,
 *      target_row) triples for Bar-Natan Gaussian elimination.
 *
 *  Key invariant (Fixed-Strand Conjecture):
 *      For fixed n strands, every vertex in G(b) has degree ≤ 6.
 *      The collapse terminates in O(poly(L)) time for fixed n.
 */

#ifndef MORSEORACLE_H
#define MORSEORACLE_H

#include <stdbool.h>

/* ================================================================
 *  Braid Representation
 * ================================================================ */

/*
 *  Braid — an n-strand braid word of length L.
 *  generators[t] ∈ {0, ..., n-2} is the Artin generator index at position t.
 *  (Corresponds to σ_{generators[t]+1} in standard notation.)
 *
 *  Note: we only handle positive braids for now. Negative crossings
 *  (σ^{-1}) can be supported by extending the graph construction.
 */
typedef struct Braid {
  int n_strands;  /* n: number of strands                     */
  int length;     /* L: number of crossings                   */
  int *generators; /* Array of size L: generator indices       */
} Braid;

Braid *Braid_create(int n_strands, int length, const int *generators);
void Braid_free(Braid *b);

/* ================================================================
 *  Augmented Rhomboid Graph
 * ================================================================ */

#define MAX_DEGREE 6 /* Proven in DegreeBound.lean */

/*
 *  AugRhomboidGraph — adjacency structure for vertices (s, t).
 *
 *  Vertices are indexed linearly: vertex(s, t) = t * n + s,
 *  where s ∈ {0, ..., n-1} and t ∈ {0, ..., L}.
 *  Total vertices: n * (L + 1).
 *
 *  For each vertex we store its neighbor list (at most MAX_DEGREE = 6).
 */
typedef struct AugRhomboidGraph {
  int n_strands;      /* n                                      */
  int n_times;        /* L + 1                                  */
  int n_vertices;     /* n * (L + 1)                            */
  int *degree;        /* Array[n_vertices]: degree of each vertex */
  int *neighbors;     /* Flat array[n_vertices * MAX_DEGREE]:
                         neighbors[v * MAX_DEGREE + k] = k-th
                         neighbor of vertex v.                   */
} AugRhomboidGraph;

/*
 *  Purpose:    Build the augmented rhomboid graph from a braid.
 *  Method:     Adds horizontal, vertical, and crossing edges
 *              per the definition in FixedStrand/Defs.lean.
 */
AugRhomboidGraph *AugRhomboidGraph_build(const Braid *b);
void AugRhomboidGraph_free(AugRhomboidGraph *g);

/*
 *  Purpose:    Convert (strand, time) to linear vertex index.
 */
static inline int vertex_index(int n_strands, int strand, int time) {
  return time * n_strands + strand;
}

/*
 *  Purpose:    Extract strand from vertex index.
 */
static inline int vertex_strand(int n_strands, int v) {
  return v % n_strands;
}

/*
 *  Purpose:    Extract time from vertex index.
 */
static inline int vertex_time(int n_strands, int v) {
  return v / n_strands;
}

/*
 *  Purpose:    Check if two vertices are adjacent in the graph.
 */
bool AugRhomboidGraph_adjacent(const AugRhomboidGraph *g, int u, int v);

/* ================================================================
 *  Independence Complex & Discrete Morse Collapse
 * ================================================================ */

/*
 *  FreePair — a free pair (σ, τ) in the independence complex.
 *  σ ⊂ τ, |τ| = |σ| + 1, and σ has no coface other than σ and τ.
 *
 *  We store both σ and τ as bitmask-encoded vertex sets.
 */
typedef struct FreePair {
  int sigma_size;     /* |σ|                                    */
  int tau_size;       /* |τ| = |σ| + 1                         */
  int *sigma_verts;   /* Sorted vertex indices of σ             */
  int *tau_verts;     /* Sorted vertex indices of τ             */
} FreePair;

/*
 *  CollapseSchedule — the deterministic sequence of free pairs
 *  found during discrete Morse collapse.
 *
 *  Each entry maps to a (chain_idx, source_col, target_row)
 *  triple for Gaussian elimination on the Khovanov complex.
 */
typedef struct CollapseSchedule {
  int count;          /* Number of collapse steps               */
  int capacity;       /* Allocated capacity                     */
  FreePair *pairs;    /* Array of free pairs                    */

  /* Khovanov mapping (filled after mapping phase): */
  int *chain_idx;     /* Which differential to reduce           */
  int *source_col;    /* Column (source generator) to eliminate */
  int *target_row;    /* Row (target generator) to eliminate    */
} CollapseSchedule;

/*
 *  Purpose:    Compute the full discrete Morse collapse schedule
 *              for the independence complex of G(b).
 *
 *  Method:     Iteratively finds free pairs and performs elementary
 *              collapses until no free pairs remain. The face set
 *              is maintained implicitly using the graph structure.
 *
 *  Complexity: O(poly(L)) for fixed n (degree ≤ 6 bounds local work).
 */
CollapseSchedule *MorseOracle_computeSchedule(const Braid *b);

void CollapseSchedule_free(CollapseSchedule *s);

/* ================================================================
 *  Full Oracle Pipeline
 * ================================================================ */

/*
 *  Purpose:    Given a braid, compute the complete elimination schedule
 *              with Khovanov chain indices filled in.
 *
 *  This is the main entry point: braid → schedule with
 *  (chain_idx, source_col, target_row) ready for blockReductionLemma.
 */
CollapseSchedule *MorseOracle_fullSchedule(const Braid *b);

#endif /* MORSEORACLE_H */
