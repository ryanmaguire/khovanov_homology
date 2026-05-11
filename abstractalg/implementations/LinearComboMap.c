#include "LinearComboMap.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "../Morphism.h"

typedef struct {
    LinearComboMap* owner;
    int index;
} LinearComboMapIteratorState;

static bool bp_is_zero(const BivariatePoly* p) {
    return p == NULL || p->num_terms == 0;
}

static BivariatePoly* clone_poly(const BivariatePoly* p) {
    if (!p) return NULL;
    BivariatePoly* copy = bp_create();
    for(int i = 0; i < p->num_terms; i++) {
        if (copy->num_terms == copy->capacity) {
            copy->capacity *= 2;
            copy->terms = realloc(copy->terms, copy->capacity * sizeof(Monomial));
        }
        copy->terms[i] = p->terms[i];
        copy->num_terms++;
    }
    return copy;
}

static bool linearComboMapIterator_hasNext(MorphismIterator* self) {
    if (!self || !self->state) return false;
    LinearComboMapIteratorState* st = self->state;
    return st->owner && st->index < st->owner->size;
}

static Morphism* linearComboMapIterator_next(MorphismIterator* self) {
    if (!self || !self->state) return NULL;
    LinearComboMapIteratorState* st = self->state;
    if (!st->owner || st->index >= st->owner->size) return NULL;
    return st->owner->terms[st->index++].morphism;
}

static void linearComboMapIterator_free(MorphismIterator* self) {
    if (!self) return;
    free(self->state);
    free(self);
}

static MorphismIterator* linearComboMapCollection_iterator(MorphismCollection* self) {
    if (!self) return NULL;
    LinearComboMap* owner = self->state;
    MorphismIterator* it = malloc(sizeof(MorphismIterator));
    LinearComboMapIteratorState* st = malloc(sizeof(LinearComboMapIteratorState));
    if (!it || !st) {
        free(it);
        free(st);
        return NULL;
    }
    st->owner = owner;
    st->index = 0;
    it->state = st;
    it->hasNext = linearComboMapIterator_hasNext;
    it->next = linearComboMapIterator_next;
    it->free = linearComboMapIterator_free;
    return it;
}

static bool linearComboMapCollection_isEmpty(MorphismCollection* self) {
    if (!self) return true;
    LinearComboMap* owner = self->state;
    return !owner || owner->size == 0;
}

static int linearComboMapCollection_size(MorphismCollection* self) {
    if (!self) return 0;
    LinearComboMap* owner = self->state;
    return owner ? owner->size : 0;
}

static MorphismCollection* LinearComboMap_terms(AbstractLinearCombo* self) {
    if (!self) return NULL;
    LinearComboMap* owner = (LinearComboMap*)self;
    MorphismCollection* col = malloc(sizeof(MorphismCollection));
    if (!col) return NULL;
    col->state = owner;
    col->iterator = linearComboMapCollection_iterator;
    col->isEmpty = linearComboMapCollection_isEmpty;
    col->size = linearComboMapCollection_size;
    return col;
}

static int morphismCompare(Morphism* a, Morphism* b) {
    if (a == b) return 0;
    if (!a) return -1;
    if (!b) return 1;
    if (a->compareTo) return a->compareTo(a, b);
    uintptr_t ua = (uintptr_t)a;
    uintptr_t ub = (uintptr_t)b;
    if (ua < ub) return -1;
    if (ua > ub) return 1;
    return 0;
}

static BivariatePoly* LinearComboMap_getCoefficient(AbstractLinearCombo* self, Morphism* term) {
    if (!self || !term) return NULL;
    LinearComboMap* owner = (LinearComboMap*)self;
    int lo = 0;
    int hi = owner->size;
    
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        int cmp = morphismCompare(owner->terms[mid].morphism, term);
        if (cmp < 0) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if (lo < owner->size && morphismCompare(owner->terms[lo].morphism, term) == 0) {
        return owner->terms[lo].coefficient;
    }
    return NULL;
}

static bool LinearComboMap_ensureCapacity(LinearComboMap* self, int minCapacity) {
    if (!self) return false;
    if (self->capacity >= minCapacity) return true;
    int newCapacity = self->capacity > 0 ? self->capacity : 1;
    while (newCapacity < minCapacity) newCapacity *= 2;
    Term* newTerms = realloc(self->terms, (size_t)newCapacity * sizeof(Term));
    if (!newTerms) return false;
    self->terms = newTerms;
    self->capacity = newCapacity;
    return true;
}

static LinearComboMap* LinearComboMap_compactImpl(LinearComboMap* self) {
    return self;
}

//purpose: initialize a new map based linear combination
//arguments: self
//argument descriptions: self is a pointer to the LinearComboMap struct to initialize
//output: void
//output description: no return value
//method: initializes the internal array and assigns function pointers
void LinearComboMap_init(LinearComboMap* self) {
    if (self) {
        self->size = 0;
        self->capacity = 10;
        self->terms = malloc(sizeof(Term) * (size_t)self->capacity);
        self->base.terms = LinearComboMap_terms;
        self->base.getCoefficient = LinearComboMap_getCoefficient;
        self->addTerm  = LinearComboMap_addTerm;
        self->addCombo = LinearComboMap_addCombo;
        self->multiply = LinearComboMap_multiply;
        self->compose  = LinearComboMap_compose;
        self->compact = LinearComboMap_compactImpl;
    }
}
//purpose: add a morphism term with an integer coefficient to the map
//arguments: self, cc, num
//argument descriptions: self is the LinearComboMap pointer, cc is the Morphism* being added, num is the int coefficient
//output: LinearComboMap*
//output description: pointer to the modified LinearComboMap
//method: we check if the term exists, update or remove it if the sum is zero, otherwise inserting the new term
LinearComboMap* LinearComboMap_addTerm(LinearComboMap* self, Morphism* cc, BivariatePoly* num) {
    if (!self || !cc || bp_is_zero(num)) {
        if (num) bp_free(num); //prevent leaks even if returning early
        return self;
    }
    int lo = 0;
    int hi = self->size;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        int cmp = morphismCompare(self->terms[mid].morphism, cc);
        if (cmp < 0) lo = mid + 1;
        else hi = mid;
    }

    //Found existing term
    if (lo < self->size && morphismCompare(self->terms[lo].morphism, cc) == 0) {
        BivariatePoly* newCoefficient = bp_add(self->terms[lo].coefficient, num);
        
        bp_free(self->terms[lo].coefficient); // Free old
        bp_free(num);                         // Free consumed argument
        
        if (bp_is_zero(newCoefficient)) {
            bp_free(newCoefficient);
            memmove(&self->terms[lo], &self->terms[lo + 1], (size_t)(self->size - lo - 1) * sizeof(Term));
            self->size--;
            if (self->compact) return self->compact(self);
            return self;
        }
        self->terms[lo].coefficient = newCoefficient;
        return self;
    }

    //Not found, insert new term
    if (!LinearComboMap_ensureCapacity(self, self->size + 1)) {
        bp_free(num); //Prevent leak
        return self; 
    }
    memmove(&self->terms[lo + 1], &self->terms[lo], (size_t)(self->size - lo) * sizeof(Term));
    self->terms[lo].morphism = cc;
    self->terms[lo].coefficient = num; 
    self->size++;
    return self;
}


//Purpose: Perform an O(M+N) linear merge on two sorted linear combinations
//Arguments: self, other
//Argument descriptions: self is the destination combo, other is the source combo to add
//Output: LinearComboMap*
//Output description: pointer to the modified, compacted LinearComboMap
//Method: Allocates a new buffer, linearly zips both sorted arrays, eliminates zeros, and swaps memory.
LinearComboMap* LinearComboMap_addCombo(LinearComboMap* self, LinearComboMap* other) {
    if (!self || !other || other->size == 0) return self;
    
    if (self->size == 0) {
        if (!LinearComboMap_ensureCapacity(self, other->size)) return self;
        for(int i = 0; i < other->size; i++) {
            self->terms[i].morphism = other->terms[i].morphism;
            self->terms[i].coefficient = clone_poly(other->terms[i].coefficient);
        }
        self->size = other->size;
        return self;
    }
    
    int maxCapacity = self->size + other->size;
    Term* mergedTerms = malloc(sizeof(Term) * maxCapacity);
    int mergedSize = 0;

    int i = 0, j = 0;
    while (i < self->size && j < other->size) {
        int cmp = morphismCompare(self->terms[i].morphism, other->terms[j].morphism);

        if (cmp == 0) {
            BivariatePoly* newCoef = bp_add(self->terms[i].coefficient, other->terms[j].coefficient);
            bp_free(self->terms[i].coefficient);
            
            if (!bp_is_zero(newCoef)) { 
                mergedTerms[mergedSize].morphism = self->terms[i].morphism;
                mergedTerms[mergedSize].coefficient = newCoef;
                mergedSize++;
            } else {
                bp_free(newCoef); 
            }
            i++; j++;
        } 
        else if (cmp < 0) {
            mergedTerms[mergedSize++] = self->terms[i++];
        } 
        else {
            mergedTerms[mergedSize].morphism = other->terms[j].morphism;
            mergedTerms[mergedSize].coefficient = clone_poly(other->terms[j].coefficient);
            mergedSize++;
            j++;
        }
    }
    while (i < self->size) {
        mergedTerms[mergedSize++] = self->terms[i++];
    }
    while (j < other->size) {
        mergedTerms[mergedSize].morphism = other->terms[j].morphism;
        mergedTerms[mergedSize].coefficient = clone_poly(other->terms[j].coefficient);
        mergedSize++; j++;
    }

    free(self->terms);
    self->terms = mergedTerms;
    self->size = mergedSize;
    self->capacity = maxCapacity;
    return self;
}
//purpose: multiply all coefficients in the combination by a scalar
//arguments: self, num
//argument descriptions: self is the LinearComboMap, num is the int scalar
//output: LinearComboMap*
//output description: pointer to the scaled LinearComboMap or null if zeroed
//method: clears the map if num is zero, otherwise iterates through terms and multiplies each coefficient
LinearComboMap* LinearComboMap_multiply(LinearComboMap* self, BivariatePoly* num) {
    if (bp_is_zero(num)) {
        if (!self) return self;
        for(int i = 0; i < self->size; i++) {
            bp_free(self->terms[i].coefficient);
        }
        self->size = 0;
        if (self->flexibleZeroLinearCombo) return self->flexibleZeroLinearCombo();
        return self;
    }
    if (!self) return self;
    for (int i = 0; i < self->size; i++) {
        BivariatePoly* newCoef = bp_multiply(self->terms[i].coefficient, num);
        bp_free(self->terms[i].coefficient);
        self->terms[i].coefficient = newCoef;
    }
    return self;
}

//purpose: vertically compose two linear combinations
//arguments: self, other
//argument descriptions: self is the first LinearComboMap, other is the LinearComboMap to compose with
//output: LinearComboMap*
//output description: a new LinearComboMap representing the composed result
//method: uses nested loops to cross multiply coefficients and compose morphisms from both combinations
LinearComboMap* LinearComboMap_compose(LinearComboMap* self, LinearComboMap* other) {
    if (!self || !other || !self->flexibleZeroLinearCombo || !self->composeMorphisms) return NULL;
    LinearComboMap* ret = self->flexibleZeroLinearCombo();
    if (!ret) return NULL;
    if (self->size == 0 || other->size == 0) return ret;
    for (int i = 0; i < self->size; i++) {
        Morphism* cc = self->terms[i].morphism;
        BivariatePoly* myCoef = self->terms[i].coefficient;
        
        for (int j = 0; j < other->size; j++) {
            Morphism* occ = other->terms[j].morphism;
            BivariatePoly* otherCoef = other->terms[j].coefficient;
            
            Morphism* composedMor = self->composeMorphisms(cc, occ);
            BivariatePoly* composedCoef = bp_multiply(myCoef, otherCoef);
            LinearComboMap_addTerm(ret, composedMor, composedCoef);
        }
    }
    if (ret && ret->compact) return ret->compact(ret);
    return ret;
}
