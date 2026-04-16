#ifndef ABSTRACT_LINEAR_COMBO_H
#define ABSTRACT_LINEAR_COMBO_H

#include <stdbool.h>
#include <stddef.h>
typedef struct Ring Ring;
typedef struct Morphism Morphism;
typedef struct MorphismIterator {
    void* state; 
    bool (*hasNext)(struct MorphismIterator* self);
    Morphism* (*next)(struct MorphismIterator* self);
    void (*free)(struct MorphismIterator* self); 
} MorphismIterator;

typedef struct MorphismCollection {
    void* state; 
    MorphismIterator* (*iterator)(struct MorphismCollection* self);
    bool (*isEmpty)(struct MorphismCollection* self);
    int (*size)(struct MorphismCollection* self);
} MorphismCollection;
typedef struct AbstractLinearCombo AbstractLinearCombo;

struct AbstractLinearCombo {
    MorphismCollection* (*terms)(AbstractLinearCombo* self);
    Ring* (*getCoefficient)(AbstractLinearCombo* self, Morphism* term);

    char* (*ringToString)(Ring* r);
    char* (*morphismToString)(Morphism* m);
};
/*
 * Purpose: Retrieves the first morphism term from the linear combination.
 * Arguments: 
 *  >self
 * Output: Morphism*
 */
Morphism* AbstractLinearCombo_firstTerm(AbstractLinearCombo* self);

/*
 * Purpose: Retrieves the coefficient of the first morphism term.
 * Arguments: 
 *  >self
 * Output: Ring*
 */
Ring* AbstractLinearCombo_firstCoefficient(AbstractLinearCombo* self);

/*
 * Purpose: Checks if the linear combination has no terms.
 * Arguments: 
 *  >self
 * Output: bool
 */
bool AbstractLinearCombo_isZero(AbstractLinearCombo* self);

/*
 * Purpose: Counts the number of terms in the linear combination.
 * Arguments: 
 *  >self
 * Output: int
 */
int AbstractLinearCombo_numberOfTerms(AbstractLinearCombo* self);

/*
 * Purpose: Generates a string representation of the linear combination.
 * Arguments: 
 *  >self
 * Output: char*
 */
char* AbstractLinearCombo_toString(AbstractLinearCombo* self);

#endif 