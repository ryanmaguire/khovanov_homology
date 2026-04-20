#include "AbstractLinearCombo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char* str;
    size_t length;
    size_t capacity;
} StringBuilder;

//Purpose: Initialize a string builder
//Arguments: sb
//Argument descriptions: Pointer to the StringBuilder struct
//Output: void
//Output description: None
//Method: Allocate 64 bytes and set the initial null terminator
static void sb_init(StringBuilder* sb) {
    sb->capacity = 64;
    sb->str = malloc(sb->capacity);
    if (sb->str) {
        sb->str[0] = '\0';
    }
    sb->length = 0;
}

//Purpose: Append text to the string builder
//Arguments: sb, text
//Argument descriptions: Pointer to builder and the string to append
//Output: void
//Output description: None
//Method: Double the memory allocation if needed and copy the new text
static void sb_append(StringBuilder* sb, const char* text) {
    if (!text || !sb->str) return;
    size_t tlen = strlen(text);
    while (sb->length + tlen + 1 > sb->capacity) {
        sb->capacity *= 2;
        sb->str = realloc(sb->str, sb->capacity);
    }
    strcpy(sb->str + sb->length, text);
    sb->length += tlen;
}

//Purpose: Retrieve the first morphism term
//Arguments: self
//Argument descriptions: Pointer to the AbstractLinearCombo instance
//Output: Morphism*
//Output description: Pointer to the first Morphism or NULL if empty
//Method: Create an iterator from terms collection, grab the first item, and free iterator
Morphism* AbstractLinearCombo_firstTerm(AbstractLinearCombo* self) {
    if (!self || !self->terms) return NULL;
    MorphismCollection* terms_col = self->terms(self);
    if (!terms_col || !terms_col->iterator) {
        if (terms_col) free(terms_col); // fix memory leak
        return NULL;
    }
    MorphismIterator* it = terms_col->iterator(terms_col);
    Morphism* first = NULL;
    if (it && it->hasNext && it->hasNext(it)) {
        first = it->next(it);
    }
    
    if (it && it->free) it->free(it); 
    free(terms_col); //fix leak
    return first;
}

//Purpose: Get the integer coefficient of the first term
//Arguments: self
//Argument descriptions: Pointer to the AbstractLinearCombo instance
//Output: int
//Output description: Integer coefficient or 0 if empty
//Method: Get the first term and pass it to the getCoefficient function
int AbstractLinearCombo_firstCoefficient(AbstractLinearCombo* self) {
    Morphism* first = AbstractLinearCombo_firstTerm(self);
    if (first != NULL && self->getCoefficient) {
        return self->getCoefficient(self, first);
    }
    return 0; 
}
//Purpose: Check if the linear combo has no terms
//Arguments: self
//Argument descriptions: Pointer to the AbstractLinearCombo instance
//Output: bool
//Output description: True if zero terms exist, false otherwise
//Method: Check the terms collection using its isEmpty function
bool AbstractLinearCombo_isZero(AbstractLinearCombo* self) {
    if (!self || !self->terms) return true;
    MorphismCollection* terms_col = self->terms(self);
    if (!terms_col || !terms_col->isEmpty) {
        if (terms_col) free(terms_col);
        return true;
    }
    
    bool result = terms_col->isEmpty(terms_col);
    free(terms_col); // fix memory leak
    return result;
}

//Purpose: Count the total number of terms
//Arguments: self
//Argument descriptions: Pointer to the AbstractLinearCombo instance
//Output: int
//Output description: Total number of topological terms
//Method: query the terms collection using its size function
int AbstractLinearCombo_numberOfTerms(AbstractLinearCombo* self) {
    if (!self || !self->terms) return 0;
    MorphismCollection* terms_col = self->terms(self);
    if (!terms_col || !terms_col->size) {
        if (terms_col) free(terms_col);
        return 0;
    }
    
    int count = terms_col->size(terms_col);
    free(terms_col); 
    return count;
}

//Purpose: Create a string of the math formula
//Arguments: self
//Argument descriptions: a pointer to the AbstractLinearCombo instance
//Output: char*
//Output description: Dynamically allocated formula string
//Method: we iterate through terms, convert coefficients, and build the string
char* AbstractLinearCombo_toString(AbstractLinearCombo* self) {
    if (!self || !self->terms || !self->morphismToString) {
        return NULL;
    }

    MorphismCollection* terms_col = self->terms(self);
    
    // Validate collection and iterator before proceeding
    if (!terms_col || !terms_col->iterator) {
        if (terms_col) free(terms_col);
        return NULL; 
    }

    MorphismIterator* i = terms_col->iterator(terms_col);
    if (!i) {
        free(terms_col);
        return NULL;
    }

    StringBuilder sb;
    sb_init(&sb);
    char coeff_buffer[32];

    if (i->hasNext && i->hasNext(i)) {
        Morphism* term = i->next(i);
        int coeff = self->getCoefficient(self, term);
        snprintf(coeff_buffer, sizeof(coeff_buffer), "%d", coeff);
        char* term_str = self->morphismToString(term);

        sb_append(&sb, coeff_buffer);
        sb_append(&sb, " * ");
        sb_append(&sb, term_str);
        free(term_str);
    }

    while (i->hasNext && i->hasNext(i)) {
        Morphism* term = i->next(i);
        int coeff = self->getCoefficient(self, term);
        snprintf(coeff_buffer, sizeof(coeff_buffer), "%d", coeff);
        char* term_str = self->morphismToString(term);

        sb_append(&sb, " + ");
        sb_append(&sb, coeff_buffer);
        sb_append(&sb, " * ");
        sb_append(&sb, term_str);
        free(term_str);
    }

    if (i->free) i->free(i);
    free(terms_col);

    return sb.str;
}