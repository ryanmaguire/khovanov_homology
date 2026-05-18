#ifndef ABSTRACT_MATRIX_H
#define ABSTRACT_MATRIX_H
#include "AbstractLinearCombo.h"

//flat array accessed via (row * num_columns) + column
typedef struct {
    int num_rows;
    int num_cols;
    AbstractLinearCombo** cells; // An array of pointers to combos
} AbstractMatrix;

//initialize the matrix with empty combos
AbstractMatrix* AbstractMatrix_create(int rows, int cols);

void AbstractMatrix_free(AbstractMatrix* self);

//get the combination at a specific cell
AbstractLinearCombo* AbstractMatrix_getEntry(AbstractMatrix* self, int row, int col);

//set the combination at a specific cell
void AbstractMatrix_putEntry(AbstractMatrix* self, int row, int col, AbstractLinearCombo* combo);

//add a new term to an existing cell
void AbstractMatrix_addTermToEntry(AbstractMatrix* self, int row, int col, Morphism* m, BivariatePoly* coeff);

#endif