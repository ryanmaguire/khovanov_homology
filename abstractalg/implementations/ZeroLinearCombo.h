#ifndef ZERO_LINEAR_COMBO_H
#define ZERO_LINEAR_COMBO_H
#include "AbstractLinearCombo.h"
typedef struct ZeroLinearCombo ZeroLinearCombo;

struct ZeroLinearCombo {
    AbstractLinearCombo base;
    ZeroLinearCombo* (*singleTermLinearCombo)(Morphism* mor, int r);
    ZeroLinearCombo* (*fixedZeroLinearCombo)(void);
    ZeroLinearCombo* (*flexibleZeroLinearCombo)(void);
};

//purpose: return the first coefficient of an empty combo
//arguments: self, a pointer to the ZeroLinearCombo instance
//output: int, since we return 0
//method: hardcoded to return 0
BivariatePoly* ZeroLinearCombo_firstCoefficient(ZeroLinearCombo* self);

//purpose: return the first morphism term
//arguments: self, a pointer to the ZeroLinearCombo instance
//output: Morphism*
//method: hardcoded to return null
Morphism* ZeroLinearCombo_firstTerm(ZeroLinearCombo* self);

//purpose: we check the total number of terms
//arguments: self, pointer to the ZeroLinearCombo instance
//output: int
//method: hardcoded to 0
int ZeroLinearCombo_numberOfTerms(ZeroLinearCombo* self);

//purpose: check if the combo is zero
//arguments: self, a pointer to the ZeroLinearCombo instance
//output: bool
//method: hardcoded to true
bool ZeroLinearCombo_isZero(ZeroLinearCombo* self);

//purpose: add another linear combo to this zero combo
//arguments: self, m. self is the ZeroLinearCombo pointer, m is the ZeroLinearCombo to add
//output: ZeroLinearCombo*, returns the passed in ZeroLinearCombo m
//method: adding to zero simply yields the other object
ZeroLinearCombo* ZeroLinearCombo_addCombo(ZeroLinearCombo* self, ZeroLinearCombo* m);

//purpose: multiply the zero combo by a scalar
//arguments: self, r, where self is the ZeroLinearCombo pointer and r is the int scalar
//output: ZeroLinearCombo*
//method: zero multiplied by anything is still zero and thus it returns itself
ZeroLinearCombo* ZeroLinearCombo_multiply(ZeroLinearCombo* self, int r);

//purpose: compose this combo with another
//arguments: self, m
//argument descriptions: self is the ZeroLinearCombo pointer, m is the other ZeroLinearCombo
//output: ZeroLinearCombo*
//output description: a new fixed zero combo
//method: composition with zero is zero so it calls 
ZeroLinearCombo* ZeroLinearCombo_compose(ZeroLinearCombo* self, ZeroLinearCombo* m);

//purpose: add a single morphism term and coefficient to the zero combo
//arguments: self, mor, r
//argument descriptions: self is the ZeroLinearCombo pointer, mor is the Morphism*, r is the int coefficient
//output: ZeroLinearCombo*
//output description: a new combo containing just the added term
ZeroLinearCombo* ZeroLinearCombo_addTerm(ZeroLinearCombo* self, Morphism* mor, int r);

//purpose: look up the coefficient of a specific term
//arguments: self, term
//argument descriptions: self is the ZeroLinearCombo pointer, term is the Morphism* to find
//output: int
//output description: always returns 0
//method: since the combo is empty the coefficient for any term is 0
BivariatePoly* ZeroLinearCombo_getCoefficient(ZeroLinearCombo* self, Morphism* term);

//purpose: retrieve the collection of terms
//arguments: self
//argument descriptions: self is the ZeroLinearCombo pointer
//output: MorphismCollection*
//output description: an empty collection object
//method: allocates and returns a collection that reports being empty
MorphismCollection* ZeroLinearCombo_terms(ZeroLinearCombo* self);

//purpose: compact the linear combination
//arguments: self
//argument descriptions: self is the ZeroLinearCombo pointer
//output: ZeroLinearCombo*
//output description: returns self
//method: a zero combo is already fully compacted so it just returns itself
ZeroLinearCombo* ZeroLinearCombo_compact(ZeroLinearCombo* self);

#endif
