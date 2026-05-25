#include "KhovanovKomplex.h"
#include "CobMatrix.h"
#include <stdlib.h>

/* Reference externally existing topology resolution functions to avoid pulling in unnecessary headers */
extern SmoothingColumn* resolution_to_smoothing_column(Knot* knot, int r);
extern CobMatrix* build_boundary_matrix(SmoothingColumn* source, SmoothingColumn* target, int r);

/* Free function used inside KnotEncoder.c */
extern void Cap_free(Cap* cap);

KhovanovKomplex* KhovanovKomplex_alloc(Knot* knot) {
    if (!knot || knot->m <= 0) {
        return NULL;
    }

    KhovanovKomplex* kh_complex = (KhovanovKomplex*)malloc(sizeof(KhovanovKomplex));
    if (!kh_complex) {
        return NULL;
    }

    kh_complex->knot = knot;
   
    /* According to current logic, total number of composite states is 2^m */
    int num_states = 1 << knot->m;
   
    /* Call Arjun's existing chain complex constructor */
    kh_complex->komplex = Komplex_create(num_states);
    if (!kh_complex->komplex) {
        free(kh_complex);
        return NULL;
    }
    return kh_complex;
}

void KhovanovKomplex_free(KhovanovKomplex* kh_complex) {
    if (!kh_complex) {
        return;
    }
   
    /* Ownership transfer: directly call Arjun's safe release interface */
    if (kh_complex->komplex) {
        Komplex_free(kh_complex->komplex);
    }
   
    /* Note: Never free kh_complex->knot. Follow the "whoever allocates, releases" principle */
    free(kh_complex);
}

bool KhovanovKomplex_build(KhovanovKomplex* kh_complex) {
    if (!kh_complex || !kh_complex->knot || !kh_complex->komplex) {
        return false;
    }

    Knot* knot = kh_complex->knot;
    Komplex* komplex = kh_complex->komplex;
    int m = knot->m;
    int num_states = 1 << m;

    /* Core assembly logic: fully reuse the loop previously validated in compute_khovanov_polynomial */
    for (int r = 0; r < num_states - 1; r++) {
       
        /* 1. Fetch topology columns and convert to your SmoothingColumn structures */
        SmoothingColumn* source = resolution_to_smoothing_column(knot, r);
        SmoothingColumn* target = resolution_to_smoothing_column(knot, r + 1);

        if (!source || !target) {
            /* Defensive memory cleanup */
            if (source) {
                for (int i = 0; i < source->n; i++) Cap_free(source->smoothings[i]);
                free(source->smoothings); free(source->numbers); free(source);
            }
            if (target) {
                for (int i = 0; i < target->n; i++) Cap_free(target->smoothings[i]);
                free(target->smoothings); free(target->numbers); free(target);
            }
            continue;
        }

        /* 2. Dispatch your existing function to generate the matrix.
           The Jordan-Wigner signs for the hypercube are already handled inside build_boundary_matrix */
        CobMatrix* boundary = build_boundary_matrix(source, target, r);
       
        /* 3. Mount the matrix directly onto Arjun's Komplex */
        komplex->differentials[r] = boundary;

        /* 4. Clean up temporary SmoothingColumn objects.
           Because build_boundary_matrix internally uses CobMatrix_create(..., false),
           it clones the matrix, so the original pointers here must be safely released. */
        for (int i = 0; i < source->n; i++) {
            if (source->smoothings[i]) Cap_free(source->smoothings[i]);
        }
        free(source->smoothings);
        free(source->numbers);
        free(source);

        for (int i = 0; i < target->n; i++) {
            if (target->smoothings[i]) Cap_free(target->smoothings[i]);
        }
        free(target->smoothings);
        free(target->numbers);
        free(target);
    }

    return true;
}
