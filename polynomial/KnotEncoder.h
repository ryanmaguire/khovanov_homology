#ifndef KNOT_ENCODER_H
#define KNOT_ENCODER_H

#include "BivariatePoly.h"
#include "IntMatrix.h"

// Simple knot representation: 4m-tuple of edge labels
typedef struct {
    int* edges;     // 4m edge labels
    int m;          // number of crossings
} Knot;

// Resolution register |r> (m bits)
typedef struct {
    int* r;         // resolution bits (0 or 1)
    int m;
} Resolution;

// Resolved knot register |K_r>
typedef struct {
    int* resolved_edges;   // reordered edge pairs after smoothing
    int m;
} ResolvedKnot;

// Loop enumeration register |L_r>
typedef struct {
    int* loop_min_labels;  // minimal edge label for each loop
    int num_loops;
} LoopRegister;

// Enhanced state |s> (labels on loops: 0=1, 1=X)
typedef struct {
    int* s;         // 0 or 1 for each loop
    int num_loops;
} EnhancedState;

// Create a hard-coded trefoil for quick testing
Knot* knot_create_trefoil(void);

// Free knot
void knot_free(Knot* k);

// Build full chain complex and compute Khovanov polynomial
BivariatePoly* compute_khovanov_polynomial(Knot* knot);

#endif