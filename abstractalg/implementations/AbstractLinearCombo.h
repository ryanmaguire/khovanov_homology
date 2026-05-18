#ifndef ABSTRACT_LINEAR_COMBO_H
#define ABSTRACT_LINEAR_COMBO_H

#include <stdbool.h>
#include <stddef.h>
#include "../../polynomial/polynomial/BivariatePoly.h"
typedef struct Morphism Morphism;
typedef struct AbstractLinearCombo AbstractLinearCombo;

typedef struct {
    Morphism* term;
    BivariatePoly* coeff;
} TermPair;

struct AbstractLinearCombo {
    int term_count;
    int capacity;
    TermPair* elements; //contiguous array of terms and coefficients
};
static inline bool AbstractLinearCombo_isZero(AbstractLinearCombo* self) {
    if (!self) return true;
    return self->term_count == 0;
}

//Purpose: Count the total number of terms
static inline int AbstractLinearCombo_numberOfTerms(AbstractLinearCombo* self) {
    if (!self) return 0;
    return self->term_count;
}

//Purpose: Retrieve the first morphism term
static inline Morphism* AbstractLinearCombo_firstTerm(AbstractLinearCombo* self) {
    if (!self || self->term_count == 0) return NULL;
    return self->elements[0].term;
}

//Purpose: Get the integer coefficient of the first term
static inline BivariatePoly* AbstractLinearCombo_firstCoefficient(AbstractLinearCombo* self) {
    if (!self || self->term_count == 0) return NULL;
    return self->elements[0].coeff;
}

//Purpose: Create a string of the math formula
char* AbstractLinearCombo_toString(AbstractLinearCombo* self, char* (*morphismToString)(Morphism*));
AbstractLinearCombo* AbstractLinearCombo_create(int initial_capacity);

//adds a term/coefficient pair, automatically resizing the array if needed
void AbstractLinearCombo_addTerm(AbstractLinearCombo* self, Morphism* term, BivariatePoly* coeff);

//frees the struct
void AbstractLinearCombo_free(AbstractLinearCombo* self);
#endif