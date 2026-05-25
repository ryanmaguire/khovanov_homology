#ifndef KHOVANOV_KOMPLEX_H
#define KHOVANOV_KOMPLEX_H

#include <stdbool.h>
#include "KnotEncoder.h"
#include "Komplex.h"

/**
 * @brief Extremely lightweight glue integration layer.
 *        Directly reuses existing Knot and Komplex pointers
 *        without introducing any extra ownership or redundant data.
 */
typedef struct {
    Knot* knot;           /* Source from the topology layer (no ownership) */
    Komplex* komplex;     /* Arjun's algebraic container (owns the data, responsible for free) */
} KhovanovKomplex;

/**
 * @brief Allocate the wrapper layer and initialize the underlying Komplex.
 * @param knot Pointer to an existing topological knot.
 */
KhovanovKomplex* KhovanovKomplex_alloc(Knot* knot);

/**
 * @brief Safely destroy the wrapper layer.
 * @note Automatically calls Komplex_free to release the algebraic layer,
 *       but strictly protects the knot from being freed.
 */
void KhovanovKomplex_free(KhovanovKomplex* kh_complex);

/**
 * @brief Core integration: uses existing parsing and boundary functions
 *        to inject topological data into Arjun's matrix chain.
 */
bool KhovanovKomplex_build(KhovanovKomplex* kh_complex);

#endif /* KHOVANOV_KOMPLEX_H */
