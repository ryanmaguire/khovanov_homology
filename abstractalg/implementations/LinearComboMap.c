#include "../Morphism.h"
#include "LinearComboMap.h"
#include <stdlib.h>
#include <string.h>

static bool bp_is_zero(const BivariatePoly* p) {
    return p == NULL || p->num_terms == 0;
}

static BivariatePoly* clone_poly(const BivariatePoly* p) {
    if (!p) return NULL;
    BivariatePoly* copy = bp_create();
    if (!copy) return NULL;
    for (int i = 0; i < p->num_terms; i++) {
        if (copy->num_terms == copy->capacity) {
            int new_capacity = copy->capacity > 0 ? copy->capacity * 2 : 1;
            Monomial* new_terms = realloc(copy->terms, (size_t)new_capacity * sizeof(Monomial));
            if (!new_terms) {
                bp_free(copy);
                return NULL;
            }
            copy->terms = new_terms;
            copy->capacity = new_capacity;
        }
        copy->terms[copy->num_terms++] = p->terms[i];
    }
    return copy;
}
LinearComboMap* LinearComboMap_create(int initial_capacity) {
    LinearComboMap* map = (LinearComboMap*)malloc(sizeof(LinearComboMap));
    if (!map) return NULL;
    
    map->combo = AbstractLinearCombo_create(initial_capacity);
    return map;
}

void LinearComboMap_free(LinearComboMap* map) {
    if (!map) return;
    if (map->combo) {
        AbstractLinearCombo_free(map->combo);
    }
    free(map);
}


void LinearComboMap_putOrAdd(LinearComboMap* map, Morphism* m, BivariatePoly* coeff) {
    if (!map || !map->combo || !m || !coeff) return;
    if (bp_is_zero(coeff)) {
        bp_free(coeff);
        return;
    }
    
    AbstractLinearCombo* combo = map->combo;
    
    for (int i = 0; i < combo->term_count; i++) {
        if (combo->elements[i].term == m) {
            
            //morphism found, so we add the polynomials together
            BivariatePoly* old_coeff = combo->elements[i].coeff;
            BivariatePoly* new_coeff = bp_add(old_coeff, coeff);
            
            //freedom
            bp_free(old_coeff);
            bp_free(coeff);

            if (bp_is_zero(new_coeff)) {
                bp_free(new_coeff);
                memmove(&combo->elements[i], &combo->elements[i + 1],
                        (size_t)(combo->term_count - i - 1) * sizeof(TermPair));
                combo->term_count--;
                return;
            }

            combo->elements[i].coeff = new_coeff;
            
            return;//we are done
        }
    }
    
    AbstractLinearCombo_addTerm(combo, m, coeff);
}

void LinearComboMap_addCombo(LinearComboMap* map, AbstractLinearCombo* other) {
    if (!map || !other) return;
    for (int i = 0; i < other->term_count; i++) {
        BivariatePoly* coeff_copy = clone_poly(other->elements[i].coeff);
        if (!coeff_copy) continue;
        LinearComboMap_putOrAdd(map, other->elements[i].term, coeff_copy);
    }
}

void LinearComboMap_multiply(LinearComboMap* self, BivariatePoly* num) {
    if (!self || !self->combo || !num) return;
    
    for (int i = 0; i < self->combo->term_count; i++) {
        BivariatePoly* old_coeff = self->combo->elements[i].coeff;
        self->combo->elements[i].coeff = bp_multiply(old_coeff, num);
        bp_free(old_coeff);
    }
}

LinearComboMap* LinearComboMap_compose(LinearComboMap* self, LinearComboMap* other, Morphism* (*composeMorphisms)(Morphism* m1, Morphism* m2)) {
    if (!self || !other || !self->combo || !other->combo || !composeMorphisms) return NULL;
    LinearComboMap* result = LinearComboMap_create(self->combo->term_count * other->combo->term_count);
    for (int i = 0; i < self->combo->term_count; i++) {
        for (int j = 0; j < other->combo->term_count; j++) {
            
            Morphism* cc = self->combo->elements[i].term;
            BivariatePoly* myCoef = self->combo->elements[i].coeff;
            
            Morphism* occ = other->combo->elements[j].term;
            BivariatePoly* otherCoef = other->combo->elements[j].coeff;
            
            Morphism* composedMor = composeMorphisms(cc, occ);
            BivariatePoly* composedCoef = bp_multiply(myCoef, otherCoef);
            
            LinearComboMap_putOrAdd(result, composedMor, composedCoef);
        }
    }
    
    return result;
}


AbstractLinearCombo* LinearComboMap_toLinearCombo(LinearComboMap* map) {
    if (!map || !map->combo) return NULL;
    
    AbstractLinearCombo* final_combo = AbstractLinearCombo_create(map->combo->term_count);
    if (!final_combo) return NULL;
    
    for (int i = 0; i < map->combo->term_count; i++) {
        BivariatePoly* coeff = map->combo->elements[i].coeff;
        
        if (!bp_is_zero(coeff)) {
            BivariatePoly* coeff_copy = clone_poly(coeff);
            if (!coeff_copy) continue;
            AbstractLinearCombo_addTerm(final_combo, map->combo->elements[i].term, coeff_copy);
        }
    }
    
    return final_combo;
}
