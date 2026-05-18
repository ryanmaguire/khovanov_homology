#ifndef LINEAR_COMBO_MAP_H
#define LINEAR_COMBO_MAP_H

#include "AbstractLinearCombo.h"
#include "../../polynomial/polynomial/BivariatePoly.h"
#include <stdbool.h>

//wrapper around array 
typedef struct {
    AbstractLinearCombo* combo; 
} LinearComboMap;

LinearComboMap* LinearComboMap_create(int initial_capacity);
void LinearComboMap_free(LinearComboMap* map);


void LinearComboMap_putOrAdd(LinearComboMap* map, Morphism* m, BivariatePoly* coeff);

void LinearComboMap_addCombo(LinearComboMap* map, AbstractLinearCombo* other);

//multiply all coefficients by a scalar polynomial
void LinearComboMap_multiply(LinearComboMap* self, BivariatePoly* num);

//cross multiply two linear combinations (used in tensor products)
LinearComboMap* LinearComboMap_compose(LinearComboMap* self, LinearComboMap* other, Morphism* (*composeMorphisms)(Morphism* m1, Morphism* m2));

//strips out any terms that cancelled to zero and returns the final
AbstractLinearCombo* LinearComboMap_toLinearCombo(LinearComboMap* map);

#endif