#ifndef ABSTRACT_LINEAR_COMBO_H
#define ABSTRACT_LINEAR_COMBO_H

#include <stdbool.h>
#include <stddef.h>

typedef struct Morphism Morphism;

typedef struct MorphismIterator {
    void* state; 
    bool (*hasNext)(struct MorphismIterator* self);
    Morphism* (*next)(struct MorphismIterator* self);
    void (*free)(struct MorphismIterator* self); 
} MorphismIterator;

typedef struct MorphismCollection {
    void* state; 
    MorphismIterator* (*iterator)(struct MorphismCollection* self);
    bool (*isEmpty)(struct MorphismCollection* self);
    int (*size)(struct MorphismCollection* self);
} MorphismCollection;

typedef struct AbstractLinearCombo AbstractLinearCombo;

struct AbstractLinearCombo {
    MorphismCollection* (*terms)(AbstractLinearCombo* self);
    int (*getCoefficient)(AbstractLinearCombo* self, Morphism* term);
    char* (*morphismToString)(Morphism* m);
};
//Purpose: Retrieve the first morphism term
//Arguments: self
//Argument descriptions: Pointer to the AbstractLinearCombo instance
//Output: Morphism*
//Output description: Pointer to the first Morphism or null if empty
//Method: Create an iterator from terms collection, grab the first item, and free iterator
Morphism* AbstractLinearCombo_firstTerm(AbstractLinearCombo* self);

//Purpose: Get the integer coefficient of the first term
//Arguments: self
//Argument descriptions: Pointer to the AbstractLinearCombo instance
//Output: int
//Output description: Integer coefficient or 0 if empty
//Method: Get the first term and pass it to the getCoefficient function
int AbstractLinearCombo_firstCoefficient(AbstractLinearCombo* self);

//Purpose: Check if the linear combo has no terms
//Arguments: self
//Argument descriptions: Pointer to the AbstractLinearCombo instance
//Output: bool
//Output description: True if zero terms exist, false otherwise
//Method: Check the terms collection using its isEmpty function
bool AbstractLinearCombo_isZero(AbstractLinearCombo* self);

//Purpose: Count the total number of terms
//Arguments: self
//Argument descriptions: Pointer to the AbstractLinearCombo instance
//Output: int
//Output description: Total number of topological terms
//Method: Query the terms collection using its size function
int AbstractLinearCombo_numberOfTerms(AbstractLinearCombo* self);

//Purpose: Create a string of the math formula
//Arguments: self
//Argument descriptions: Pointer to the AbstractLinearCombo instance
//Output: char*
// Output description: Dynamically allocated formula string
//Method: Iterate through terms, convert coefficients, and build the string dynamically
char* AbstractLinearCombo_toString(AbstractLinearCombo* self);

#endif