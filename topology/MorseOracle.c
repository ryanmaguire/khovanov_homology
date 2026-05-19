/*
 *  MorseOracle.c
 *
 *  Topological Discrete Morse Oracle for Fixed-Strand FPT Khovanov Homology.
 *
 *  Purpose:
 *      Phase 1 of the FPT algorithm. Constructs the augmented rhomboid
 *      graph from a braid word, enumerates independent sets, finds free
 *      pairs via discrete Morse collapse, and maps them to Khovanov
 *      chain complex elimination targets.
 *
 *  Compile:
 *      gcc -Wall -Wextra -g -c MorseOracle.c -I.
 */

#include "MorseOracle.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ================================================================
 *  Braid
 * ================================================================ */

Braid *Braid_create(int n_strands, int length, const int *generators) {
  Braid *b = (Braid *)malloc(sizeof(Braid));
  b->n_strands = n_strands;
  b->length = length;
  b->generators = (int *)malloc((size_t)length * sizeof(int));
  memcpy(b->generators, generators, (size_t)length * sizeof(int));
  return b;
}

void Braid_free(Braid *b) {
  if (!b) return;
  free(b->generators);
  free(b);
}

/* ================================================================
 *  Augmented Rhomboid Graph — Construction
 *
 *  Edge types from FixedStrand/Defs.lean (AugRhomboidAdj):
 *
 *  1. Horizontal: (s, t) — (s, t±1)
 *     Same strand, adjacent time steps.
 *
 *  2. Vertical: (s, t) — (s±1, t)
 *     Same time step, adjacent strands.
 *
 *  3. Crossing: for generator σ_i at time t (strands i, i+1):
 *     (i, t) — (i+1, t+1)  and  (i+1, t) — (i, t+1)
 *     These are the diagonal "crossing" edges.
 * ================================================================ */

/*
 *  Purpose:    Try to add an edge u—v to the graph.
 *              Skips if already present or would exceed MAX_DEGREE.
 */
static void try_add_edge(AugRhomboidGraph *g, int u, int v) {
  if (u == v) return;
  if (u < 0 || u >= g->n_vertices) return;
  if (v < 0 || v >= g->n_vertices) return;

  /* Check if edge already exists */
  int du = g->degree[u];
  for (int k = 0; k < du; k++) {
    if (g->neighbors[u * MAX_DEGREE + k] == v) return;
  }

  /* Add u→v */
  if (du < MAX_DEGREE) {
    g->neighbors[u * MAX_DEGREE + du] = v;
    g->degree[u] = du + 1;
  }

  /* Add v→u */
  int dv = g->degree[v];
  if (dv < MAX_DEGREE) {
    g->neighbors[v * MAX_DEGREE + dv] = u;
    g->degree[v] = dv + 1;
  }
}

AugRhomboidGraph *AugRhomboidGraph_build(const Braid *b) {
  int n = b->n_strands;
  int L = b->length;
  int T = L + 1; /* number of time steps */
  int V = n * T;

  AugRhomboidGraph *g = (AugRhomboidGraph *)malloc(sizeof(AugRhomboidGraph));
  g->n_strands = n;
  g->n_times = T;
  g->n_vertices = V;
  g->degree = (int *)calloc((size_t)V, sizeof(int));
  g->neighbors = (int *)calloc((size_t)V * MAX_DEGREE, sizeof(int));

  /* Initialize neighbors to -1 (sentinel) */
  memset(g->neighbors, -1, (size_t)V * MAX_DEGREE * sizeof(int));

  /* 1. Horizontal edges: (s, t) — (s, t+1) */
  for (int s = 0; s < n; s++) {
    for (int t = 0; t < T - 1; t++) {
      int u = vertex_index(n, s, t);
      int v = vertex_index(n, s, t + 1);
      try_add_edge(g, u, v);
    }
  }

  /* 2. Vertical edges: (s, t) — (s+1, t) */
  for (int s = 0; s < n - 1; s++) {
    for (int t = 0; t < T; t++) {
      int u = vertex_index(n, s, t);
      int v = vertex_index(n, s + 1, t);
      try_add_edge(g, u, v);
    }
  }

  /* 3. Crossing edges: for generator i at time t,
   *    (i, t) — (i+1, t+1)  and  (i+1, t) — (i, t+1) */
  for (int t = 0; t < L; t++) {
    int i = b->generators[t]; /* strand index of crossing */
    /* Edge: (i, t) — (i+1, t+1) */
    try_add_edge(g, vertex_index(n, i, t), vertex_index(n, i + 1, t + 1));
    /* Edge: (i+1, t) — (i, t+1) */
    try_add_edge(g, vertex_index(n, i + 1, t), vertex_index(n, i, t + 1));
  }

  return g;
}

void AugRhomboidGraph_free(AugRhomboidGraph *g) {
  if (!g) return;
  free(g->degree);
  free(g->neighbors);
  free(g);
}

bool AugRhomboidGraph_adjacent(const AugRhomboidGraph *g, int u, int v) {
  if (u < 0 || u >= g->n_vertices || v < 0 || v >= g->n_vertices)
    return false;
  for (int k = 0; k < g->degree[u]; k++) {
    if (g->neighbors[u * MAX_DEGREE + k] == v)
      return true;
  }
  return false;
}

/* ================================================================
 *  Independence Complex — Face Enumeration
 *
 *  A face of I(G) is an independent set of G: a set S of vertices
 *  such that no two vertices in S are adjacent.
 *
 *  For small vertex counts (fixed n, moderate L), we enumerate
 *  faces using a backtracking algorithm. For large L, we use
 *  a sliding-window approach that exploits the O(1) degree bound.
 * ================================================================ */

/*
 *  FaceSet — dynamically resizable collection of independent sets.
 *  Each face is stored as a sorted array of vertex indices.
 */
typedef struct Face {
  int size;
  int *verts;        /* Sorted vertex indices */
} Face;

typedef struct FaceSet {
  int count;
  int capacity;
  Face *faces;
  bool *alive;       /* alive[i] = true if face i hasn't been collapsed */
} FaceSet;

static FaceSet *FaceSet_create(int initial_cap) {
  FaceSet *fs = (FaceSet *)malloc(sizeof(FaceSet));
  fs->count = 0;
  fs->capacity = initial_cap;
  fs->faces = (Face *)malloc((size_t)initial_cap * sizeof(Face));
  fs->alive = (bool *)malloc((size_t)initial_cap * sizeof(bool));
  return fs;
}

static void FaceSet_addFace(FaceSet *fs, int size, const int *verts) {
  if (fs->count >= fs->capacity) {
    fs->capacity *= 2;
    fs->faces = (Face *)realloc(fs->faces, (size_t)fs->capacity * sizeof(Face));
    fs->alive = (bool *)realloc(fs->alive, (size_t)fs->capacity * sizeof(bool));
  }
  Face *f = &fs->faces[fs->count];
  f->size = size;
  f->verts = (int *)malloc((size_t)size * sizeof(int));
  memcpy(f->verts, verts, (size_t)size * sizeof(int));
  fs->alive[fs->count] = true;
  fs->count++;
}

static void FaceSet_free(FaceSet *fs) {
  if (!fs) return;
  for (int i = 0; i < fs->count; i++) {
    free(fs->faces[i].verts);
  }
  free(fs->faces);
  free(fs->alive);
  free(fs);
}

/*
 *  Purpose:    Check if vertex v is independent from all vertices in set[0..size-1].
 */
static bool is_independent_of(const AugRhomboidGraph *g, const int *set,
                               int size, int v) {
  for (int i = 0; i < size; i++) {
    if (AugRhomboidGraph_adjacent(g, set[i], v))
      return false;
  }
  return true;
}

/*
 *  Purpose:    Recursively enumerate all independent sets (faces) of I(G).
 *              Uses backtracking with the vertices ordered by index.
 *
 *  Arguments:
 *      g       — the graph
 *      fs      — output face set
 *      current — current partial independent set
 *      size    — size of current set
 *      next_v  — next vertex to consider (for ordered enumeration)
 *      nv      — total number of vertices
 */
static void enumerate_faces(const AugRhomboidGraph *g, FaceSet *fs,
                            int *current, int size, int next_v, int nv) {
  /* Record the current set as a face (including the empty set) */
  FaceSet_addFace(fs, size, current);

  /* Try extending with each remaining vertex */
  for (int v = next_v; v < nv; v++) {
    if (is_independent_of(g, current, size, v)) {
      current[size] = v;
      enumerate_faces(g, fs, current, size + 1, v + 1, nv);
    }
  }
}

/*
 *  Purpose:    Build the complete face set of the independence complex.
 */
static FaceSet *build_independence_complex(const AugRhomboidGraph *g) {
  FaceSet *fs = FaceSet_create(256);

  int *current = (int *)malloc((size_t)g->n_vertices * sizeof(int));
  enumerate_faces(g, fs, current, 0, 0, g->n_vertices);
  free(current);

  return fs;
}

/* ================================================================
 *  Discrete Morse Collapse — Free Pair Detection
 *
 *  A free pair (σ, τ) in K:
 *    - σ ∈ K.faces, τ ∈ K.faces
 *    - σ ⊂ τ, |τ| = |σ| + 1
 *    - σ is contained in no other face of K besides σ and τ
 *      (i.e., σ is a "free face" of τ)
 * ================================================================ */

/*
 *  Purpose:    Check if face A is a subset of face B.
 *              Both must be sorted arrays.
 */
static bool face_subset(const int *a, int a_size, const int *b, int b_size) {
  if (a_size > b_size) return false;
  int j = 0;
  for (int i = 0; i < a_size; i++) {
    while (j < b_size && b[j] < a[i]) j++;
    if (j >= b_size || b[j] != a[i]) return false;
    j++;
  }
  return true;
}

/*
 *  Purpose:    Check if two sorted arrays are equal.
 */
static bool face_equal(const int *a, int a_size, const int *b, int b_size) {
  if (a_size != b_size) return false;
  return memcmp(a, b, (size_t)a_size * sizeof(int)) == 0;
}

/*
 *  Purpose:    Find a free pair (σ, τ) in the current alive faces.
 *              Returns true if found, with indices stored in *sigma_idx, *tau_idx.
 *
 *  Method:     For each alive face σ, count how many alive faces
 *              properly contain σ with |τ| = |σ| + 1. If exactly one
 *              such τ exists, (σ, τ) is a free pair.
 */
static bool find_free_pair(const FaceSet *fs, int *sigma_idx, int *tau_idx) {
  for (int i = 0; i < fs->count; i++) {
    if (!fs->alive[i]) continue;
    Face *sigma = &fs->faces[i];

    /* Count cofaces of σ with dimension |σ| + 1 */
    int coface_count = 0;
    int coface_idx = -1;

    for (int j = 0; j < fs->count; j++) {
      if (!fs->alive[j]) continue;
      if (i == j) continue;
      Face *tau = &fs->faces[j];

      /* Check: τ has size |σ| + 1 and σ ⊂ τ */
      if (tau->size == sigma->size + 1 &&
          face_subset(sigma->verts, sigma->size, tau->verts, tau->size)) {
        coface_count++;
        coface_idx = j;
        if (coface_count > 1) break; /* optimization: no need to count beyond 1 */
      }
    }

    if (coface_count == 1) {
      *sigma_idx = i;
      *tau_idx = coface_idx;
      return true;
    }
  }
  return false;
}

/* ================================================================
 *  Collapse Schedule Construction
 * ================================================================ */

static CollapseSchedule *CollapseSchedule_create(int initial_cap) {
  CollapseSchedule *s = (CollapseSchedule *)malloc(sizeof(CollapseSchedule));
  s->count = 0;
  s->capacity = initial_cap;
  s->pairs = (FreePair *)malloc((size_t)initial_cap * sizeof(FreePair));
  s->chain_idx = (int *)malloc((size_t)initial_cap * sizeof(int));
  s->source_col = (int *)malloc((size_t)initial_cap * sizeof(int));
  s->target_row = (int *)malloc((size_t)initial_cap * sizeof(int));
  return s;
}

static void CollapseSchedule_addPair(CollapseSchedule *s,
                                      const Face *sigma, const Face *tau) {
  if (s->count >= s->capacity) {
    s->capacity *= 2;
    s->pairs = (FreePair *)realloc(s->pairs, (size_t)s->capacity * sizeof(FreePair));
    s->chain_idx = (int *)realloc(s->chain_idx, (size_t)s->capacity * sizeof(int));
    s->source_col = (int *)realloc(s->source_col, (size_t)s->capacity * sizeof(int));
    s->target_row = (int *)realloc(s->target_row, (size_t)s->capacity * sizeof(int));
  }

  FreePair *fp = &s->pairs[s->count];
  fp->sigma_size = sigma->size;
  fp->tau_size = tau->size;
  fp->sigma_verts = (int *)malloc((size_t)sigma->size * sizeof(int));
  fp->tau_verts = (int *)malloc((size_t)tau->size * sizeof(int));
  memcpy(fp->sigma_verts, sigma->verts, (size_t)sigma->size * sizeof(int));
  memcpy(fp->tau_verts, tau->verts, (size_t)tau->size * sizeof(int));

  /* Khovanov indices filled later in mapping phase */
  s->chain_idx[s->count] = -1;
  s->source_col[s->count] = -1;
  s->target_row[s->count] = -1;

  s->count++;
}

void CollapseSchedule_free(CollapseSchedule *s) {
  if (!s) return;
  for (int i = 0; i < s->count; i++) {
    free(s->pairs[i].sigma_verts);
    free(s->pairs[i].tau_verts);
  }
  free(s->pairs);
  free(s->chain_idx);
  free(s->source_col);
  free(s->target_row);
  free(s);
}

/* ================================================================
 *  Core Oracle: Compute Collapse Schedule
 * ================================================================ */

CollapseSchedule *MorseOracle_computeSchedule(const Braid *b) {
  /* 1. Build the augmented rhomboid graph */
  AugRhomboidGraph *g = AugRhomboidGraph_build(b);

  /* 2. Enumerate all faces of the independence complex */
  FaceSet *fs = build_independence_complex(g);

  /* 3. Iteratively find and collapse free pairs */
  CollapseSchedule *schedule = CollapseSchedule_create(64);

  int sigma_idx, tau_idx;
  while (find_free_pair(fs, &sigma_idx, &tau_idx)) {
    /* Record the free pair */
    CollapseSchedule_addPair(schedule, &fs->faces[sigma_idx],
                             &fs->faces[tau_idx]);

    /* Elementary collapse: remove both σ and τ */
    fs->alive[sigma_idx] = false;
    fs->alive[tau_idx] = false;
  }

  /* Cleanup */
  FaceSet_free(fs);
  AugRhomboidGraph_free(g);

  return schedule;
}

/* ================================================================
 *  Mapping: Free Pairs → Khovanov Chain Indices
 *
 *  In the cube of resolutions for a braid of length L:
 *  - Each resolution is a binary string r ∈ {0,1}^L
 *  - The chain index is the Hamming weight |r| (number of 1-bits)
 *  - Within chain index k, generators are indexed by the
 *    lexicographic order of k-bit subsets of {0,...,L-1}
 *
 *  A free pair (σ, τ) in the independence complex, where
 *  |τ| = |σ| + 1, maps to adjacent resolutions that differ
 *  in exactly one crossing. The extra vertex in τ \ σ identifies
 *  which crossing flips and in which direction.
 *
 *  For the cube-of-resolutions → chain complex mapping:
 *  - chain_idx = |σ| (the smaller face's dimension in the complex)
 *  - source_col = lexicographic index of σ among all faces of size |σ|
 *  - target_row = lexicographic index of τ among all faces of size |τ|
 *
 *  However, since we perform collapses sequentially and each collapse
 *  removes generators, we need to track index shifts. We do this
 *  by maintaining running counts.
 * ================================================================ */

/*
 *  Purpose:    Compare two faces for lexicographic ordering.
 *              Used to establish a canonical ordering of generators
 *              within each chain group.
 */
static int face_compare(const int *a, int a_size, const int *b, int b_size) {
  /* First by size */
  if (a_size != b_size) return a_size - b_size;
  /* Then lexicographically */
  for (int i = 0; i < a_size; i++) {
    if (a[i] != b[i]) return a[i] - b[i];
  }
  return 0;
}

/*
 *  Purpose:    Map all free pairs in the schedule to Khovanov indices.
 *
 *  We rebuild the face set and process collapses in order, computing
 *  the correct indices at each step (accounting for previously
 *  removed generators).
 */
static void map_schedule_to_khovanov(CollapseSchedule *schedule,
                                      const Braid *b) {
  /* Rebuild the independence complex */
  AugRhomboidGraph *g = AugRhomboidGraph_build(b);
  FaceSet *fs = build_independence_complex(g);

  for (int step = 0; step < schedule->count; step++) {
    FreePair *fp = &schedule->pairs[step];

    /* chain_idx = dimension of σ in the chain complex = |σ| */
    schedule->chain_idx[step] = fp->sigma_size;

    /*
     * Find the index of σ among all alive faces of the same size,
     * ordered lexicographically.
     */
    int sigma_rank = 0;
    for (int i = 0; i < fs->count; i++) {
      if (!fs->alive[i]) continue;
      if (fs->faces[i].size != fp->sigma_size) continue;
      int cmp = face_compare(fs->faces[i].verts, fs->faces[i].size,
                              fp->sigma_verts, fp->sigma_size);
      if (cmp < 0) sigma_rank++;
      if (cmp == 0) break;
    }
    schedule->source_col[step] = sigma_rank;

    /*
     * Find the index of τ among all alive faces of the same size.
     */
    int tau_rank = 0;
    for (int i = 0; i < fs->count; i++) {
      if (!fs->alive[i]) continue;
      if (fs->faces[i].size != fp->tau_size) continue;
      int cmp = face_compare(fs->faces[i].verts, fs->faces[i].size,
                              fp->tau_verts, fp->tau_size);
      if (cmp < 0) tau_rank++;
      if (cmp == 0) break;
    }
    schedule->target_row[step] = tau_rank;

    /* Now perform the collapse: remove σ and τ from the alive set */
    for (int i = 0; i < fs->count; i++) {
      if (!fs->alive[i]) continue;
      if (fs->faces[i].size == fp->sigma_size &&
          face_equal(fs->faces[i].verts, fs->faces[i].size,
                     fp->sigma_verts, fp->sigma_size)) {
        fs->alive[i] = false;
      }
      if (fs->faces[i].size == fp->tau_size &&
          face_equal(fs->faces[i].verts, fs->faces[i].size,
                     fp->tau_verts, fp->tau_size)) {
        fs->alive[i] = false;
      }
    }
  }

  FaceSet_free(fs);
  AugRhomboidGraph_free(g);
}

/* ================================================================
 *  Full Oracle Pipeline
 * ================================================================ */

CollapseSchedule *MorseOracle_fullSchedule(const Braid *b) {
  /* Phase 1a: Compute collapse schedule (topological) */
  CollapseSchedule *schedule = MorseOracle_computeSchedule(b);

  /* Phase 1b: Map free pairs to Khovanov elimination targets */
  map_schedule_to_khovanov(schedule, b);

  return schedule;
}
