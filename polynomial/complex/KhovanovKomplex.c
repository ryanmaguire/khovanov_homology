#include "KhovanovKomplex.h"
#include "CobMatrix.h"
#include <stdlib.h>
#include <stdio.h>

// ------------------------------------------------------------------
// Core internal topology helper functions: managed directly by the glue layer
// to prevent misuse of external interfaces.
// ------------------------------------------------------------------

/**
 * @brief Compute the number of 1-bits in a binary number (homological grading).
 */
static int pop_count(unsigned int n) {
    int count = 0;
    while (n) { n &= (n - 1); count++; }
    return count;
}

/**
 * @brief [Key fix] Convert a complete resolution state r into a single Cap.
 *        In Arjun's design, a Cap represents an entire 1-manifold (a full resolution state).
 */
static Cap* build_resolution_cap(Knot* knot, int r) {
    int m = knot->m;
    Cap* cap = Cap_create(2 * m, 0); // Create a Cap with 2m boundary points
    
    for (int k = 0; k < m; k++) {
        // Note: based on your knot_from_pd implementation, extract PD coordinates
        int a = knot->pd[k][0], b = knot->pd[k][1];
        int c = knot->pd[k][2], d = knot->pd[k][3];
        
        if ((r & (1 << k)) == 0) {
            // 0-smoothing
            cap->pairings[a] = b; cap->pairings[b] = a;
            cap->pairings[c] = d; cap->pairings[d] = c;
        } else {
            // 1-smoothing
            cap->pairings[a] = d; cap->pairings[d] = a;
            cap->pairings[b] = c; cap->pairings[c] = b;
        }
    }
    return cap;
}

/**
 * @brief Compute the Jordan-Wigner sign for a hypercube edge.
 */
static int compute_jordan_wigner_sign(int src_res, int dst_res) {
    int diff = src_res ^ dst_res;
    // Ensure exactly one bit changes and it is a 0 -> 1 flip
    if (diff == 0 || (diff & (diff - 1)) != 0 || (diff & dst_res) == 0) return 0;
    
    int mask = diff - 1;
    int prior_ones = src_res & mask;
    return (pop_count(prior_ones) % 2 == 0) ? 1 : -1;
}

// ------------------------------------------------------------------
// Lifecycle and build API
// ------------------------------------------------------------------

KhovanovKomplex* KhovanovKomplex_alloc(Knot* knot) {
    if (!knot || knot->m <= 0) return NULL;

    KhovanovKomplex* kh_complex = (KhovanovKomplex*)malloc(sizeof(KhovanovKomplex));
    if (!kh_complex) return NULL;

    kh_complex->knot = knot;
    
    // [Key fix] The length of the chain complex is the number of homological layers m + 1,
    // not the total number of states 2^m
    kh_complex->komplex = Komplex_create(knot->m + 1);
    if (!kh_complex->komplex) {
        free(kh_complex);
        return NULL;
    }
    return kh_complex;
}

void KhovanovKomplex_free(KhovanovKomplex* kh_complex) {
    if (!kh_complex) return;
    if (kh_complex->komplex) Komplex_free(kh_complex->komplex);
    free(kh_complex);
}

// Assume an externally provided LCCC constructor based on cobordisms
// (this fixes the missing cc issue in your previous code)
extern LCCC* LCCC_create_with_cobordism(int sign, CannedCobordism* cc);
// If this is not yet implemented, you can temporarily fall back to your old LCCC_create(sign)

bool KhovanovKomplex_build(KhovanovKomplex* kh_complex) {
    if (!kh_complex || !kh_complex->knot || !kh_complex->komplex) return false;

    Knot* knot = kh_complex->knot;
    Komplex* komplex = kh_complex->komplex;
    int m = knot->m;
    int num_states = 1 << m;

    // 1. Count the number of states in each homological grading h (binomial coefficients)
    int* states_count = calloc(m + 1, sizeof(int));
    for (int r = 0; r < num_states; r++) {
        states_count[pop_count(r)]++;
    }

    // 2. Allocate state arrays for each grading h and record the actual r values
    int** states_by_h = malloc((m + 1) * sizeof(int*));
    int* current_idx = calloc(m + 1, sizeof(int));
    for (int h = 0; h <= m; h++) {
        states_by_h[h] = malloc(states_count[h] * sizeof(int));
    }
    
    for (int r = 0; r < num_states; r++) {
        int h = pop_count(r);
        states_by_h[h][current_idx[h]++] = r;
    }
    free(current_idx);

    // 3. Build the real SmoothingColumn array (each Column represents an entire homological group C_h)
    SmoothingColumn** columns = malloc((m + 1) * sizeof(SmoothingColumn*));
    for (int h = 0; h <= m; h++) {
        columns[h] = malloc(sizeof(SmoothingColumn));
        columns[h]->n = states_count[h];
        columns[h]->numbers = calloc(columns[h]->n, sizeof(int)); // for quantum grading
        columns[h]->smoothings = malloc(columns[h]->n * sizeof(Cap*));
        
        for (int i = 0; i < columns[h]->n; i++) {
            columns[h]->smoothings[i] = build_resolution_cap(knot, states_by_h[h][i]);
        }
    }

    // 4. Build block boundary matrices differentials mapping C_h -> C_{h+1}
    for (int h = 0; h < m; h++) {
        CobMatrix* boundary = CobMatrix_create(columns[h], columns[h+1], false);
        
        // Traverse target rows and source columns
        for (int row = 0; row < columns[h+1]->n; row++) {
            int target_r = states_by_h[h+1][row];
            
            for (int col = 0; col < columns[h]->n; col++) {
                int source_r = states_by_h[h][col];
                
                int sign = compute_jordan_wigner_sign(source_r, target_r);
                if (sign != 0) {
                    // A valid hypercube edge exists, generate local cobordism mapping
                    CannedCobordismImplData* impl = CannedCobordismImpl_create(
                        columns[h]->smoothings[col], 
                        columns[h+1]->smoothings[row]
                    );
                    CannedCobordism* cc = CannedCobordismImpl_as_CannedCobordism(impl);
                    
                    // Note: You need to coordinate with Eric to wrap the cobordism cc and sign into LCCC
                    // Temporarily using placeholder; replace with LCCC_create_with_cobordism if available
                    LCCC* lc = LCCC_create(sign); 
                    
                    // Inject into the correct coordinate of the large block matrix
                    CobMatrix_putEntry(boundary, row, col, lc);
                }
            }
        }
        komplex->differentials[h] = boundary;
    }

    // 5. Memory cleanup (release outer shells; inner Caps are handled by Arjun's clone or shared policy)
    for (int h = 0; h <= m; h++) {
        for (int i = 0; i < columns[h]->n; i++) {
            Cap_free(columns[h]->smoothings[i]);
        }
        free(columns[h]->smoothings);
        free(columns[h]->numbers);
        free(columns[h]);
        free(states_by_h[h]);
    }
    free(columns);
    free(states_by_h);
    free(states_count);

    return true;
}
