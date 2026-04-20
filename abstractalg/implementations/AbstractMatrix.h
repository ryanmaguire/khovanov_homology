#ifndef ABSTRACT_MATRIX_H
#define ABSTRACT_MATRIX_H

#include <stdbool.h>
#include <stddef.h>

typedef struct Morphism Morphism;

typedef struct AbstractMatrix AbstractMatrix;

struct AbstractMatrix {
    Morphism* (*getEntry)(AbstractMatrix* self, int row, int column);
    void (*putEntry)(AbstractMatrix* self, int row, int column, Morphism* t);
    Morphism* (*addMorphism)(Morphism* m1, Morphism* m2);
};

//Purpose: safely add a morphism to a specific matrix cell
//Arguments: self, row, column, t
//Argument descriptions: self is a pointer to the AbstractMatrix instance, row and column are the integer coordinates, t is the Morphism* to insert
//Output: void
//Output description: no return value
//Method: retrieve the existing cell entry, if it is NULL insert t directly, otherwise combine them using addMorphism and insert the result
void AbstractMatrix_addEntry(AbstractMatrix* self, int row, int column, Morphism* t);

#endif