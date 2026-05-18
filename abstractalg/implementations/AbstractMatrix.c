#include "AbstractMatrix.h"
#include "LinearComboMap.h"
#include <stdlib.h>

AbstractMatrix* AbstractMatrix_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return NULL;
    AbstractMatrix* m = (AbstractMatrix*)malloc(sizeof(AbstractMatrix));
    if (!m) return NULL;
    m->num_rows = rows;
    m->num_cols = cols;
    //the flat array of pointers
    m->cells = (AbstractLinearCombo**)malloc(rows * cols * sizeof(AbstractLinearCombo*));
    if (!m->cells) {
        free(m);
        return NULL;
    }
    //initialize every cell with an empty linear combination
    for (int i = 0; i < rows * cols; i++) {
        m->cells[i] = AbstractLinearCombo_create(0);
    }
    return m;
}

void AbstractMatrix_free(AbstractMatrix* self) {
    if (!self) return;
    if (self->cells) {
        for (int i = 0; i < self->num_rows * self->num_cols; i++) {
            if (self->cells[i]) {
                AbstractLinearCombo_free(self->cells[i]);
            }
        }
        free(self->cells);
    }
    free(self); //free
}

AbstractLinearCombo* AbstractMatrix_getEntry(AbstractMatrix* self, int row, int col) {
    if (!self || row < 0 || row >= self->num_rows || col < 0 || col >= self->num_cols) return NULL;
    
    int index = (row * self->num_cols) + col;
    return self->cells[index];
}

void AbstractMatrix_putEntry(AbstractMatrix* self, int row, int col, AbstractLinearCombo* combo) {
    if (!self || row < 0 || row >= self->num_rows || col < 0 || col >= self->num_cols) return;
    
    int index = (row * self->num_cols) + col;
    if (self->cells[index] && self->cells[index] != combo) {
        AbstractLinearCombo_free(self->cells[index]);
    }
    
    self->cells[index] = combo;
}

void AbstractMatrix_addTermToEntry(AbstractMatrix* self, int row, int col, Morphism* m, BivariatePoly* coeff) {
    if (!self || row < 0 || row >= self->num_rows || col < 0 || col >= self->num_cols) return;
    
    int index = (row * self->num_cols) + col;
    
    AbstractLinearCombo* cellCombo = self->cells[index];
    LinearComboMap mapWrapper;
    mapWrapper.combo = cellCombo;
    
    LinearComboMap_putOrAdd(&mapWrapper, m, coeff);
}