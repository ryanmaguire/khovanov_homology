#include "ZeroLinearCombo.h"
#include <stdlib.h>

//purpose: provide hasNext logic for an empty iterator
//arguments: self
//argument descriptions: self is a pointer to the MorphismIterator
//output: bool
//output description: always returns false
//method: hardcoded to false
static bool emptyIterator_hasNext(MorphismIterator* self) {
    return false;
}

//purpose: provide next logic for an empty iterator
//arguments: self
//argument descriptions: self is a pointer to the MorphismIterator
//output: Morphism*
//output description: always returns null
//method: hardcoded to null
static Morphism* emptyIterator_next(MorphismIterator* self) {
    return NULL;
}

//purpose: clean up memory for the empty iterator
//arguments: self
//argument descriptions: self is a pointer to the MorphismIterator
//output: void
//output description: no return value
//method: calls standard free
static void emptyIterator_free(MorphismIterator* self) {
    if (self) free(self);
}

//purpose: spawn an iterator for an empty collection
//arguments: self
//argument descriptions: self is a pointer to the MorphismCollection
//output: MorphismIterator*
//output description: a freshly allocated iterator wired to do nothing
//method: allocates memory and wires up the empty logic function pointers
static MorphismIterator* emptyCollection_iterator(MorphismCollection* self) {
    MorphismIterator* it = malloc(sizeof(MorphismIterator));
    if (it) {
        it->state = NULL;
        it->hasNext = emptyIterator_hasNext;
        it->next = emptyIterator_next;
        it->free = emptyIterator_free;
    }
    return it;
}

//purpose: provide isEmpty logic for the empty collection
//arguments: self
//argument descriptions: self is a pointer to the MorphismCollection
//output: bool
//output description: always true
//method: hardcoded to true
static bool emptyCollection_isEmpty(MorphismCollection* self) {
    return true;
}

//purpose: provide size logic for the empty collection
//arguments: self
//argument descriptions: self is a pointer to the MorphismCollection
//output: int
//output description: always 0
//method: hardcoded to 0
static int emptyCollection_size(MorphismCollection* self) {
    return 0;
}

static BivariatePoly* ZeroLinearCombo_zeroPoly(void) {
    static Monomial terms[1];
    static BivariatePoly poly = {
        .terms = terms,
        .num_terms = 0,
        .capacity = 1,
    };
    poly.num_terms = 0;
    return &poly;
}

BivariatePoly* ZeroLinearCombo_firstCoefficient(ZeroLinearCombo* self) {
    return ZeroLinearCombo_zeroPoly();
}

Morphism* ZeroLinearCombo_firstTerm(ZeroLinearCombo* self) {
    return NULL;
}

int ZeroLinearCombo_numberOfTerms(ZeroLinearCombo* self) {
    return 0;
}

bool ZeroLinearCombo_isZero(ZeroLinearCombo* self) {
    return true;
}

ZeroLinearCombo* ZeroLinearCombo_addCombo(ZeroLinearCombo* self, ZeroLinearCombo* m) {
    return m;
}

ZeroLinearCombo* ZeroLinearCombo_multiply(ZeroLinearCombo* self, int r) {
    return self;
}

ZeroLinearCombo* ZeroLinearCombo_compose(ZeroLinearCombo* self, ZeroLinearCombo* m) {
    if (self && self->fixedZeroLinearCombo) {
        return self->fixedZeroLinearCombo();
    }
    return NULL;
}

ZeroLinearCombo* ZeroLinearCombo_addTerm(ZeroLinearCombo* self, Morphism* mor, int r) {
    if (self && self->singleTermLinearCombo) {
        return self->singleTermLinearCombo(mor, r);
    }
    return NULL;
}

BivariatePoly* ZeroLinearCombo_getCoefficient(ZeroLinearCombo* self, Morphism* term) {
    return ZeroLinearCombo_zeroPoly();
}

MorphismCollection* ZeroLinearCombo_terms(ZeroLinearCombo* self) {
    MorphismCollection* col = malloc(sizeof(MorphismCollection));
    if (col) {
        col->state = NULL;
        col->iterator = emptyCollection_iterator;
        col->isEmpty = emptyCollection_isEmpty;
        col->size = emptyCollection_size;
    }
    return col;
}

ZeroLinearCombo* ZeroLinearCombo_compact(ZeroLinearCombo* self) {
    return self;
}
