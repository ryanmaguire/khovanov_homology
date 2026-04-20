#include "LinearComboMap.h"
#include <stdlib.h>

//purpose: initialize a new map based linear combination
//arguments: self
//argument descriptions: self is a pointer to the LinearComboMap struct to initialize
//output: void
//output description: no return value
//method: allocates or sets up the internal coefficients map and assigns function pointers
void LinearComboMap_init(LinearComboMap* self) {
    //ensure validity
    if (self) {
        self->coefficients = NULL; 
    }
}

//purpose: add a morphism term with an integer coefficient to the map
//arguments: self, cc, num
//argument descriptions: self is the LinearComboMap pointer, cc is the Morphism* being added, num is the int coefficient
//output: LinearComboMap*
//output description: pointer to the modified LinearComboMap
//method: we check if the term exists, update or remove it if the sum is zero, otherwise inserting the new term
LinearComboMap* LinearComboMap_addTerm(LinearComboMap* self, Morphism* cc, int num) {
    if (!self || !self->coefficients) return self;

    MorphismIntMap* map = self->coefficients;
    //preexisting or not
    if (map->containsKey(map, cc)) {
        int newCoefficient = map->get(map, cc) + num;
        //remove if cancels out
        if (newCoefficient== 0) {
            map->remove(map, cc);
            if (self->compact) return self->compact(self);
        } else{
            //update with new coefficient
            map->put(map, cc, newCoefficient);
        }
    } else {
        if (num != 0) {
            map->put(map, cc, num);
        }
    }
    return self;
}

//purpose: add the contents of another linear combination to this one
//arguments: self, other
//argument descriptions: self is the destination LinearComboMap, other is the source LinearComboMap to add
//output: LinearComboMap*
//output description: pointer to the compacted LinearComboMap after addition
//method: we iterate through the terms of other and call addTerm for each one
LinearComboMap* LinearComboMap_addCombo(LinearComboMap* self, LinearComboMap* other) {
    //validity check
    if (!self || !other || !other->coefficients || !other->coefficients->keys) return self;
    MorphismCollection* terms = other->coefficients->keys(other->coefficients);
    if (!terms) return self;
    if (terms->iterator) {
        MorphismIterator* it = terms->iterator(terms);
        while (it && it->hasNext && it->hasNext(it)) {
            Morphism* cc = it->next(it);
            int coef = other->coefficients->get(other->coefficients, cc);
            //add into linear combination
            LinearComboMap_addTerm(self, cc, coef);
        }
        //clean up
        if (it && it->free) it->free(it);
    }
    free(terms);
    if (self->compact) return self->compact(self);
    return self;
}

//purpose: multiply all coefficients in the combination by a scalar
//arguments: self, num
//argument descriptions: self is the LinearComboMap, num is the int scalar
//output: LinearComboMap*
//output description: pointer to the scaled LinearComboMap or null if zeroed
//method: clears the map if num is zero, otherwise iterates through terms and multiplies each coefficient
LinearComboMap* LinearComboMap_multiply(LinearComboMap* self, int num) {
    if (!self || !self->coefficients) return self;

    MorphismIntMap* map = self->coefficients;

    if (num == 0) {
        //multiplying by zero clears
        if (map->clear) map->clear(map);
        if (self->flexibleZeroLinearCombo) return self->flexibleZeroLinearCombo();
        return self; 
    } else {
        if (!map->keys) return self;
        MorphismCollection* terms = map->keys(map);
        if (!terms) return self;
        if (terms->iterator) {
            MorphismIterator* it = terms->iterator(terms);
            //loop and scale
            while (it && it->hasNext && it->hasNext(it)) {
                Morphism* cc = it->next(it);
                int currentCoef = map->get(map, cc);
                map->put(map, cc, currentCoef * num);
            }
            if (it && it->free) it->free(it);
        }
        free(terms);
        return self;
    }
}

//purpose: vertically compose two linear combinations
//arguments: self, other
//argument descriptions: self is the first LinearComboMap, other is the LinearComboMap to compose with
//output: LinearComboMap*
//output description: a new LinearComboMap representing the composed result
//method: uses nested loops to cross multiply coefficients and compose morphisms from both combinations
LinearComboMap* LinearComboMap_compose(LinearComboMap* self, LinearComboMap* other) {
    if (!self || !other || !self->flexibleZeroLinearCombo || !self->composeMorphisms) return NULL;
    if (!self->coefficients || !self->coefficients->keys || !other->coefficients || !other->coefficients->keys) return NULL;
    MorphismCollection* myTerms = self->coefficients->keys(self->coefficients);
    MorphismCollection* otherTerms = other->coefficients->keys(other->coefficients);
    if (!myTerms || !otherTerms) {
        if (myTerms) free(myTerms);
        if (otherTerms) free(otherTerms);
        return self->flexibleZeroLinearCombo();
    }
    bool selfZero = myTerms->isEmpty ? myTerms->isEmpty(myTerms) : true;
    bool otherZero = otherTerms->isEmpty ? otherTerms->isEmpty(otherTerms) : true;
    if (selfZero || otherZero) {
        free(myTerms);
        free(otherTerms);
        return self->flexibleZeroLinearCombo();
    }
    LinearComboMap* ret = self->flexibleZeroLinearCombo();
    if (myTerms->iterator && otherTerms->iterator) {
        MorphismIterator* myIt = myTerms->iterator(myTerms);
        while (myIt && myIt->hasNext && myIt->hasNext(myIt)) {
            Morphism* cc = myIt->next(myIt);
            int myCoef = self->coefficients->get(self->coefficients, cc);
            MorphismIterator* otherIt = otherTerms->iterator(otherTerms);
            while (otherIt && otherIt->hasNext && otherIt->hasNext(otherIt)) {
                Morphism* occ = otherIt->next(otherIt);
                int otherCoef = other->coefficients->get(other->coefficients, occ);    
                Morphism* composedMor = self->composeMorphisms(cc, occ);
                //crossmultiply
                int composedCoef = myCoef * otherCoef;
                LinearComboMap_addTerm(ret, composedMor, composedCoef);
            }
            if (otherIt && otherIt->free) otherIt->free(otherIt);
        }
        if (myIt && myIt->free) myIt->free(myIt);
    }
    free(myTerms);
    free(otherTerms);
    if (ret && ret->compact) return ret->compact(ret);
    return ret;
}