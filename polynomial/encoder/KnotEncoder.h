#ifndef KNOT_ENCODER_H
#define KNOT_ENCODER_H

#include "../polynomial/BivariatePoly.h"

//knot representation: 4m-tuple of edge labels
typedef struct {
    int** pd;       //planar Diagram codes (a 2D array of crossings)
    int writhe;     //the writhe of the knot
    int* edges;     //4m edge labels
    int m;          //number of crossings
} Knot;

//resolution register |r> (m bits)
typedef struct {
    int* r;         //resolution bits (0 or 1)
    int m;
} Resolution;

//resolved knot register |K_r>
typedef struct {
    int* resolved_edges;   //reordered edge pairs after smoothing
    int m;
} ResolvedKnot;

//loop enumeration register |L_r>
typedef struct {
    int* loop_min_labels;  //minimal edge label for each loop
    int num_loops;
} LoopRegister;

//enhanced state |s> (labels on loops: 0=1, 1=X)
typedef struct {
    int* s;         //0 or 1 for each loop
    int num_loops;
} EnhancedState;

//we hard code a trefoil for quick testing
Knot* knot_create_trefoil(void);
//helper to create knot from PD code
Knot* knot_from_pd(int pd[][4], int m);

int knot_count_loops(Knot* knot, int r);

//free knot
void knot_free(Knot* k);

//build full chain complex and compute Khovanov polynomial
BivariatePoly* compute_khovanov_polynomial(Knot* knot);

#endif
