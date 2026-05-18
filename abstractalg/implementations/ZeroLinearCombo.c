#include "ZeroLinearCombo.h"
#include <stdlib.h>

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

ZeroLinearCombo* ZeroLinearCombo_create(void) {
    ZeroLinearCombo* z = malloc(sizeof(ZeroLinearCombo));
    if (!z) return NULL;
    z->base.term_count = 0;
    z->base.capacity = 0;
    z->base.elements = NULL;
    return z;
}

void ZeroLinearCombo_free(ZeroLinearCombo* self) {
    free(self);
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
