#include "AbstractMatrix.h"

//Purpose: add a morphism to a specific matrix cell
//Arguments: self, row, column, t
//Argument descriptions: self is a pointer to the AbstractMatrix instance, row and column are the integer coordinates, t is the Morphism* to insert
//Output: void
//Output description: no return value
//Method: retrieve the existing cell entry, if it is null we insert t directly, otherwise we combine them using addMorphism and insert the result
void AbstractMatrix_addEntry(AbstractMatrix* self, int row, int column, Morphism* t) {
    if (!self || !self->getEntry || !self->putEntry || !self->addMorphism) {
        return;
    }

    Morphism* existingEntry = self->getEntry(self, row, column);
    
    if (existingEntry == NULL) {
        self->putEntry(self, row, column, t);
    } else {
        Morphism* combined = self->addMorphism(existingEntry, t);
        self->putEntry(self, row, column, combined);
    }
}