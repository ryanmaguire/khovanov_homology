#include "KhovanovKomplex.h"
#include "../../topology/CobMatrix.h"
#include "../../topology/CannedCobordismImpl.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Macro to safely isolate structural access. Protects against future changes in Arjun's Komplex. */
#define SET_COMPLEX_DIFFERENTIAL(k_ptr, idx, mat_ptr) ((k_ptr)->differentials[(idx)] = (mat_ptr))

//debug helper
static void debug_report_khbuild(const char *location, const char *msg, int h, int row, int col) {
    char command[2048];
    snprintf(command, sizeof(command),
             "sh -lc 'source .dbg/trefoil-poincare-hang.env 2>/dev/null || { DEBUG_SERVER_URL=\"http://127.0.0.1:7777/event\"; DEBUG_SESSION_ID=\"trefoil-poincare-hang\"; }; curl -sX POST \"$DEBUG_SERVER_URL\" -H \"Content-Type: application/json\" -d \"{\\\"sessionId\\\":\\\"$DEBUG_SESSION_ID\\\",\\\"runId\\\":\\\"pre-fix\\\",\\\"hypothesisId\\\":\\\"B\\\",\\\"location\\\":\\\"%s\\\",\\\"msg\\\":\\\"%s\\\",\\\"data\\\":{\\\"h\\\":%d,\\\"row\\\":%d,\\\"col\\\":%d}}\" >/dev/null 2>&1'",
             location, msg, h, row, col);
    (void)system(command);
}
//debug end

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
    int cap_size = 4 * m;
    Cap* cap = Cap_create(cap_size, 0); 
    if (!cap) return NULL; 

    int max_edge = 0;
    for (int k = 0; k < m; k++) {
        for (int i = 0; i < 4; i++) {
            if (knot->pd[k][i] > max_edge) max_edge = knot->pd[k][i];
        }
    }
    
    int* first_port = (int*)malloc((max_edge + 1) * sizeof(int));
    for (int i = 0; i <= max_edge; i++) first_port[i] = -1;
    
    int* external_mate = (int*)malloc(cap_size * sizeof(int));
    
    for (int k = 0; k < m; k++) {
        for (int i = 0; i < 4; i++) {
            int port = 4 * k + i;
            int edge = knot->pd[k][i];
            if (first_port[edge] == -1) {
                first_port[edge] = port;
            } else {
                int mate = first_port[edge];
                external_mate[port] = mate;
                external_mate[mate] = port;
            }
        }
    }
    free(first_port);

    for (int k = 0; k < m; k++) {
        int p0 = 4 * k + 0, p1 = 4 * k + 1;
        int p2 = 4 * k + 2, p3 = 4 * k + 3;
        
        int in_mate[4];
        if ((r & (1 << k)) == 0) {
            in_mate[0] = p1; in_mate[1] = p0;
            in_mate[2] = p3; in_mate[3] = p2;
        } else {
            in_mate[0] = p3; in_mate[3] = p0;
            in_mate[1] = p2; in_mate[2] = p1;
        }
        
        cap->pairings[p0] = external_mate[in_mate[0]];
        cap->pairings[p1] = external_mate[in_mate[1]];
        cap->pairings[p2] = external_mate[in_mate[2]];
        cap->pairings[p3] = external_mate[in_mate[3]];
    }
    
    free(external_mate);
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
    kh_complex->states_by_h = (int**)calloc(m + 1, sizeof(int*)); 
    kh_complex->columns = (SmoothingColumn**)calloc(m + 1, sizeof(SmoothingColumn*));

    if (!kh_complex->states_count || !kh_complex->states_by_h || !kh_complex->columns) {
        goto alloc_failure;
    }

    //count how many states belong to each homological degre
    for (int r = 0; r < num_states; r++) {
        kh_complex->states_count[pop_count(r)]++;
    }

    
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

    //precomputation
    for (int h = 0; h <= m; h++) {
        kh_complex->columns[h] = (SmoothingColumn*)malloc(sizeof(SmoothingColumn));
        if (!kh_complex->columns[h]) goto alloc_failure;

        kh_complex->columns[h]->n = kh_complex->states_count[h];
        kh_complex->columns[h]->numbers = (int*)calloc(kh_complex->columns[h]->n, sizeof(int)); 
        kh_complex->columns[h]->smoothings = (Cap**)malloc(kh_complex->columns[h]->n * sizeof(Cap*));
        
        if (!kh_complex->columns[h]->numbers || !kh_complex->columns[h]->smoothings) goto alloc_failure;

        for (int i = 0; i < kh_complex->columns[h]->n; i++) {
            //build cap 
            int r = kh_complex->states_by_h[h][i];
            kh_complex->columns[h]->smoothings[i] = build_resolution_cap(knot, r);
            kh_complex->columns[h]->numbers[i] = knot_count_loops(knot, r); 
            
            if (!kh_complex->columns[h]->smoothings[i]) goto alloc_failure;
        }
    }

    return kh_complex;

alloc_failure:
    //free if malloc doesnt work
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
    //debug start
    debug_report_khbuild("KhovanovKomplex.c:build:enter", "start", m, 0, 0);
    //debug end

    /* Iterate through adjacent homological degrees (h -> h+1) to build block matrices */
    for (int h = 0; h < m; h++) {
        //debug block
        debug_report_khbuild("KhovanovKomplex.c:block:start", "block", h, kh_complex->columns[h+1]->n, kh_complex->columns[h]->n);
        //debug end
        
        CobMatrix* boundary = CobMatrix_create(kh_complex->columns[h], kh_complex->columns[h+1], false);
        if (!boundary) return false;
        
        /* Double loop mapping every source state in C_h to every target state in C_{h+1} */
        for (int row = 0; row < kh_complex->columns[h+1]->n; row++) {
            int target_r = kh_complex->states_by_h[h+1][row];
            if (row == 0 || row == kh_complex->columns[h+1]->n - 1) {
                //debug row
                debug_report_khbuild("KhovanovKomplex.c:block:row", "row", h, row, target_r);
                //debug end
            }
            
            for (int col = 0; col < kh_complex->columns[h]->n; col++) {
                int source_r = kh_complex->states_by_h[h][col];
                
                /* Calculate hypercube edge parity */
                int sign = compute_jordan_wigner_sign(source_r, target_r);
                
                /* If a valid topological edge exists (sign != 0) */
                if (sign != 0) {
                    if (h == 0) {
                        //debug edge
                        debug_report_khbuild("KhovanovKomplex.c:h0-edge:before-create", "edge", h, row, col);
                        //debug end
                    }
                    if (h == 0 && row == 0) {
                        //debug first
                        debug_report_khbuild("KhovanovKomplex.c:first-edge:before-create", "first", h, row, col);
                        //debug end
                    }
                    /* Create the geometric surface (cobordism) spanning the two states */
                    CannedCobordismImplData* impl = CannedCobordismImpl_create(
                        kh_complex->columns[h]->smoothings[col], 
                        kh_complex->columns[h+1]->smoothings[row]
                    );
                    if (h == 0) {
                        //debug made
                        debug_report_khbuild("KhovanovKomplex.c:h0-edge:after-create", impl ? "made" : "null", h, row, col);
                        //debug end
                    }
                    if (h == 0 && row == 0) {
                        //debug first done
                        debug_report_khbuild("KhovanovKomplex.c:first-edge:after-create", impl ? "first ok" : "first null", h, row, col);
                        //debug end
                    }
                    CannedCobordism* cc = CannedCobordismImpl_as_CannedCobordism(impl);
                    
                    /* Wrap the geometry and the algebraic sign into a linear combination object */
                    LCCC* lc = LCCC_create_with_cobordism(sign, cc); 
                    if (h == 0 && row == 0) {
                        //debug lccc
                        debug_report_khbuild("KhovanovKomplex.c:first-edge:after-lccc", lc ? "lccc ok" : "lccc null", h, row, col);
                        //debug end
                    }
                    
                    /* Inject into the algebraic grid */
                    CobMatrix_putEntry(boundary, row, col, lc);
                    if (h == 0 && row == 0) {
                        //debug put
                        debug_report_khbuild("KhovanovKomplex.c:first-edge:after-put", "put", h, row, col);
                        //debug end
                    }
                }
            }
        }
        
        /* Mount the constructed block matrix into the chain complex */
        SET_COMPLEX_DIFFERENTIAL(komplex, h, boundary);
        //debug done
        debug_report_khbuild("KhovanovKomplex.c:block:done", "done", h, kh_complex->columns[h+1]->n, kh_complex->columns[h]->n);
        //debug end
    }

    //debug exit
    debug_report_khbuild("KhovanovKomplex.c:build:exit", "exit", m, 0, 0);
    //debug end
    return true;
}
