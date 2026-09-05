#ifndef KHOVANOV_KOMPLEX_H
#define KHOVANOV_KOMPLEX_H

#include <stdbool.h>
#include "../encoder/KnotEncoder.h"
#include "../../topology/Komplex.h"

/**
 * @struct KhovanovKomplex
 * @brief Extremely lightweight, highly optimized glue layer.
 * It bridges the topological knot resolutions 
 * and the algebraic chain complex 
 */
typedef struct {
    Knot* knot;                   /* Pointer to the topology source (No ownership, DO NOT FREE) */
    Komplex* komplex;           
    
    /* --- Performance Optimization: Pre-computed & Cached Topology Data --- */
    /* By computing these during alloc(), the build() phase becomes incredibly fast and avoids memory fragmentation. */
    
    int* states_count;            /* states_count[h] stores the total number of valid resolutions at homological degree h. */
    int** states_by_h;            /* states_by_h[h][i] stores the actual resolution bitmask (r) for the i-th state at degree h. */
    SmoothingColumn** columns;    /* columns[h] caches the entire SmoothingColumn (all Caps) for homological degree h. */

} KhovanovKomplex;

/**
 * @brief Allocates the wrapper layer, pre-computes combinatorics, and caches all topology columns.
 * @param knot Pointer to an existing topological knot.
 * @return A fully initialized KhovanovKomplex pointer, or NULL if memory allocation fails.
 */
KhovanovKomplex* KhovanovKomplex_alloc(Knot* knot);

/**
 * @brief Safely destroys the wrapper layer and all cached topological templates.
 * @note Automatically calls Komplex_free to release the algebraic layer, but strictly protects the 'knot' pointer.
 * @param kh_complex The wrapper instance to free.
 */
void KhovanovKomplex_free(KhovanovKomplex* kh_complex);

/**
 * @param kh_complex The initialized wrapper instance.
 * @return true on successful build, false otherwise.
 */
bool KhovanovKomplex_build(KhovanovKomplex* kh_complex);

#endif /* KHOVANOV_KOMPLEX_H */
