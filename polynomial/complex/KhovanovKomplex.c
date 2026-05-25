#include "KhovanovKomplex.h"
#include "CobMatrix.h"
#include "CannedCobordismImpl.h" /* Required for CannedCobordismImpl_create and casting */
#include <stdlib.h>
#include <stdio.h>

/* Macro to safely isolate structural access. Protects against future changes in Arjun's Komplex. */
#define SET_COMPLEX_DIFFERENTIAL(k_ptr, idx, mat_ptr) ((k_ptr)->differentials[(idx)] = (mat_ptr))

/* External definitions provided by Arjun/Eric's algebraic layer */
extern LCCC* LCCC_create_with_cobordism(int sign, CannedCobordism* cc);
extern void Cap_free(Cap* cap);

/* ========================================================================= */
/* 1. INTERNAL TOPOLOGICAL HELPER FUNCTIONS                                  */
/* ========================================================================= */

/**
 * @brief Computes the number of set bits (1s) in an integer.
 * In Khovanov Homology, the number of 1-smoothings determines the homological degree (h).
 */
static int pop_count(unsigned int n) {
    int count = 0;
    while (n) { 
        n &= (n - 1); 
        count++; 
    }
    return count;
}

/**
 * @brief Converts a complete resolution bitmask 'r' into a single topological Cap.
 */
static Cap* build_resolution_cap(Knot* knot, int r) {
    int m = knot->m;
    Cap* cap = Cap_create(2 * m, 0); 
    if (!cap) return NULL; 
    
    for (int k = 0; k < m; k++) {
        int a = knot->pd[k][0], b = knot->pd[k][1];
        int c = knot->pd[k][2], d = knot->pd[k][3];
        
        if ((r & (1 << k)) == 0) { 
            /* 0-smoothing: connect (a,b) and (c,d) */
            cap->pairings[a] = b; cap->pairings[b] = a;
            cap->pairings[c] = d; cap->pairings[d] = c;
        } else { 
            /* 1-smoothing: connect (a,d) and (b,c) */
            cap->pairings[a] = d; cap->pairings[d] = a;
            cap->pairings[b] = c; cap->pairings[c] = b;
        }
    }
    return cap;
}

/**
 * @brief Computes the Jordan-Wigner sign for the boundary operator on the hypercube.
 * @return +1 or -1 if it's a valid edge (exactly one 0->1 flip), otherwise 0 (no connection).
 */
static int compute_jordan_wigner_sign(int src_res, int dst_res) {
    int diff = src_res ^ dst_res;
    
    /* Ensure exactly one bit changed, and it changed from 0 to 1 (source -> target) */
    if (diff == 0 || (diff & (diff - 1)) != 0 || (diff & dst_res) == 0) {
        return 0;
    }
    
    /* Count the number of 1s appearing BEFORE the flipped bit to determine parity */
    int mask = diff - 1;
    int prior_ones = src_res & mask;
    return (pop_count(prior_ones) % 2 == 0) ? 1 : -1;
}

/* ========================================================================= */
/* 2. LIFECYCLE API (ALLOC & FREE)                                           */
/* ========================================================================= */

KhovanovKomplex* KhovanovKomplex_alloc(Knot* knot) {
    if (!knot || knot->m <= 0) return NULL;

    KhovanovKomplex* kh_complex = (KhovanovKomplex*)malloc(sizeof(KhovanovKomplex));
    if (!kh_complex) return NULL;

    int m = knot->m;
    int num_states = 1 << m;

    /* Initialize pointers to NULL to allow safe partial cleanup upon failure */
    kh_complex->knot = knot;
    kh_complex->komplex = NULL;
    kh_complex->states_count = NULL;
    kh_complex->states_by_h = NULL;
    kh_complex->columns = NULL;

    /* Create the algebraic container. Chain length is m + 1 (degrees 0 to m). */
    kh_complex->komplex = Komplex_create(m + 1);
    if (!kh_complex->komplex) goto alloc_failure;

    /* --- Pre-compute Combinatorics --- */
    kh_complex->states_count = (int*)calloc(m + 1, sizeof(int));
    kh_complex->states_by_h = (int**)calloc(m + 1, sizeof(int*)); /* Crucial: use calloc for safe freeing */
    kh_complex->columns = (SmoothingColumn**)calloc(m + 1, sizeof(SmoothingColumn*));

    if (!kh_complex->states_count || !kh_complex->states_by_h || !kh_complex->columns) {
        goto alloc_failure;
    }

    kh_complex->columns[h]->numbers[i] = knot_count_loops(knot, r);

    /* Count how many states belong to each homological degree */
    for (int r = 0; r < num_states; r++) {
        kh_complex->states_count[pop_count(r)]++;
    }

    /* Allocate sub-arrays and distribute the resolution bitmasks (r) into their respective degree buckets */
    int* current_idx = (int*)calloc(m + 1, sizeof(int));
    if (!current_idx) goto alloc_failure;

    for (int h = 0; h <= m; h++) {
        if (kh_complex->states_count[h] > 0) {
            kh_complex->states_by_h[h] = (int*)malloc(kh_complex->states_count[h] * sizeof(int));
            if (!kh_complex->states_by_h[h]) {
                free(current_idx);
                goto alloc_failure;
            }
        }
    }
    
    for (int r = 0; r < num_states; r++) {
        int h = pop_count(r);
        kh_complex->states_by_h[h][current_idx[h]++] = r;
    }
    free(current_idx);

    /* --- Pre-compute all SmoothingColumns (Topology Caching) --- */
    for (int h = 0; h <= m; h++) {
        kh_complex->columns[h] = (SmoothingColumn*)malloc(sizeof(SmoothingColumn));
        if (!kh_complex->columns[h]) goto alloc_failure;

        kh_complex->columns[h]->n = kh_complex->states_count[h];
        kh_complex->columns[h]->numbers = (int*)calloc(kh_complex->columns[h]->n, sizeof(int)); 
        kh_complex->columns[h]->smoothings = (Cap**)malloc(kh_complex->columns[h]->n * sizeof(Cap*));
        
        if (!kh_complex->columns[h]->numbers || !kh_complex->columns[h]->smoothings) goto alloc_failure;

        for (int i = 0; i < kh_complex->columns[h]->n; i++) {
            kh_complex->columns[h]->smoothings[i] = build_resolution_cap(knot, kh_complex->states_by_h[h][i]);
            if (!kh_complex->columns[h]->smoothings[i]) goto alloc_failure;
        }
    }

    return kh_complex;

alloc_failure:
    /* Elegant degradation: cleanly free everything if any malloc fails */
    KhovanovKomplex_free(kh_complex);
    return NULL;
}

void KhovanovKomplex_free(KhovanovKomplex* kh_complex) {
    if (!kh_complex) return;
    
    /* 1. Free algebraic layer */
    if (kh_complex->komplex) {
        Komplex_free(kh_complex->komplex);
    }
    
    /* Variables needed to navigate the cleanup */
    int m = (kh_complex->knot != NULL) ? kh_complex->knot->m : 0;
    
    /* 2. Free cached topological columns and their internal Caps */
    if (kh_complex->columns) {
        for (int h = 0; h <= m; h++) {
            if (kh_complex->columns[h]) {
                if (kh_complex->columns[h]->smoothings) {
                    for (int i = 0; i < kh_complex->columns[h]->n; i++) {
                        if (kh_complex->columns[h]->smoothings[i]) {
                            Cap_free(kh_complex->columns[h]->smoothings[i]);
                        }
                    }
                    free(kh_complex->columns[h]->smoothings);
                }
                if (kh_complex->columns[h]->numbers) {
                    free(kh_complex->columns[h]->numbers);
                }
                free(kh_complex->columns[h]);
            }
        }
        free(kh_complex->columns);
    }
    
    /* 3. Free combinatorics caching */
    if (kh_complex->states_by_h) {
        for (int h = 0; h <= m; h++) {
            if (kh_complex->states_by_h[h]) free(kh_complex->states_by_h[h]);
        }
        free(kh_complex->states_by_h);
    }
    
    if (kh_complex->states_count) {
        free(kh_complex->states_count);
    }
    
    /* 4. Free the wrapper itself. We DO NOT free kh_complex->knot. */
    free(kh_complex);
}

/* ========================================================================= */
/* 3. CORE BUILD ENGINE                                                      */
/* ========================================================================= */

bool KhovanovKomplex_build(KhovanovKomplex* kh_complex) {
    if (!kh_complex || !kh_complex->knot || !kh_complex->komplex || !kh_complex->columns) {
        return false;
    }

    Komplex* komplex = kh_complex->komplex;
    int m = kh_complex->knot->m;

    /* Iterate through adjacent homological degrees (h -> h+1) to build block matrices */
    for (int h = 0; h < m; h++) {
        
        CobMatrix* boundary = CobMatrix_create(kh_complex->columns[h], kh_complex->columns[h+1], false);
        if (!boundary) return false;
        
        /* Double loop mapping every source state in C_h to every target state in C_{h+1} */
        for (int row = 0; row < kh_complex->columns[h+1]->n; row++) {
            int target_r = kh_complex->states_by_h[h+1][row];
            
            for (int col = 0; col < kh_complex->columns[h]->n; col++) {
                int source_r = kh_complex->states_by_h[h][col];
                
                /* Calculate hypercube edge parity */
                int sign = compute_jordan_wigner_sign(source_r, target_r);
                
                /* If a valid topological edge exists (sign != 0) */
                if (sign != 0) {
                    /* Create the geometric surface (cobordism) spanning the two states */
                    CannedCobordismImplData* impl = CannedCobordismImpl_create(
                        kh_complex->columns[h]->smoothings[col], 
                        kh_complex->columns[h+1]->smoothings[row]
                    );
                    CannedCobordism* cc = CannedCobordismImpl_as_CannedCobordism(impl);
                    
                    /* Wrap the geometry and the algebraic sign into a linear combination object */
                    LCCC* lc = LCCC_create_with_cobordism(sign, cc); 
                    
                    /* Inject into the algebraic grid */
                    CobMatrix_putEntry(boundary, row, col, lc);
                }
            }
        }
        
        /* Mount the constructed block matrix into the chain complex */
        SET_COMPLEX_DIFFERENTIAL(komplex, h, boundary);
    }

    return true;
}
