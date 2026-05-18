#ifndef ZERO_LINEAR_COMBO_H
#define ZERO_LINEAR_COMBO_H
#include "AbstractLinearCombo.h"
typedef struct ZeroLinearCombo ZeroLinearCombo;

struct ZeroLinearCombo {
    AbstractLinearCombo base;
};

ZeroLinearCombo* ZeroLinearCombo_create(void);
void ZeroLinearCombo_free(ZeroLinearCombo* self);

BivariatePoly* ZeroLinearCombo_firstCoefficient(ZeroLinearCombo* self);

Morphism* ZeroLinearCombo_firstTerm(ZeroLinearCombo* self);

int ZeroLinearCombo_numberOfTerms(ZeroLinearCombo* self);

bool ZeroLinearCombo_isZero(ZeroLinearCombo* self);

#endif
