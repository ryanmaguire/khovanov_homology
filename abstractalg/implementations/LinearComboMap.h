#ifndef LINEAR_COMBO_MAP_H
#define LINEAR_COMBO_MAP_H

#include "AbstractLinearCombo.h"
#include "../../polynomial/polynomial/BivariatePoly.h"
#include <stdbool.h>
typedef struct {
    Morphism* morphism;
    BivariatePoly* coefficient; //Now a polynomial pointer
} Term;

typedef struct LinearComboMap LinearComboMap;

struct LinearComboMap {
    AbstractLinearCombo base;
    Term* terms;
    int size;
    int capacity;
    LinearComboMap* (*addTerm)(LinearComboMap* self, Morphism* cc, BivariatePoly* num);
    LinearComboMap* (*addCombo)(LinearComboMap* self, LinearComboMap* other);
    LinearComboMap* (*multiply)(LinearComboMap* self, BivariatePoly* num);
    LinearComboMap* (*compose)(LinearComboMap* self, LinearComboMap* other);
    LinearComboMap* (*compact)(LinearComboMap* self);
    LinearComboMap* (*flexibleZeroLinearCombo)(void);
    Morphism* (*composeMorphisms)(Morphism* m1, Morphism* m2);
};

//purpose: initialize a new map based linear combination
//arguments: self
//argument descriptions: self is a pointer to the LinearComboMap struct to initialize
//output: void
//output description: no return value
//method: initializes the internal array and assigns function pointers
void LinearComboMap_init(LinearComboMap* self);

//purpose: add a morphism term with an integer coefficient to the map
//arguments: self, cc, num
//argument descriptions: self is the LinearComboMap pointer, cc is the Morphism* being added, num is the int coefficient
//output: LinearComboMap*, a pointer to the modified LinearComboMap
//method: checks if the term exists, updates or removes it if the sum is zero, otherwise inserts the new term
LinearComboMap* LinearComboMap_addTerm(LinearComboMap* self, Morphism* cc, BivariatePoly* num);

//purpose: add the contents of another linear combination to this one
//arguments: self, other
//argument descriptions: self is the destination LinearComboMap, other is the source LinearComboMap to add
//output: LinearComboMap*
//output description: pointer to the compacted LinearComboMap after addition
//method: iterates through the terms of other and calls addTerm for each one
LinearComboMap* LinearComboMap_addCombo(LinearComboMap* self, LinearComboMap* other);

//purpose: multiply all coefficients in the combination by a scalar
//arguments: self, num
//argument descriptions: self is the LinearComboMap, num is the int scalar
//output: LinearComboMap*
//output description: pointer to the scaled LinearComboMap or NULL if zeroed
//method: clears the map if num is zero, otherwise iterates through terms and multiplies each coefficient
LinearComboMap* LinearComboMap_multiply(LinearComboMap* self, BivariatePoly* num);

//purpose: vertically compose two linear combinations
//arguments: self, other
//argument descriptions: self is the first LinearComboMap, other is the LinearComboMap to compose with
//output: LinearComboMap*
//output description: a new LinearComboMap representing the composed result
//method: uses nested loops to cross multiply coefficients and compose morphisms from both combinations
LinearComboMap* LinearComboMap_compose(LinearComboMap* self, LinearComboMap* other);

#endif
