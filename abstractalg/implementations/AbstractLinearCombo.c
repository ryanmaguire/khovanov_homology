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

//Purpose: Create a string of the math formula
//Arguments: self, morphismToString
//Argument descriptions: a pointer to the AbstractLinearCombo, and a function pointer to stringify a Morphism
//Output: char*
//Output description: Dynamically allocated formula string
//Method: Iterate through the flat elements array, convert coefficients, and build the string
char* AbstractLinearCombo_toString(AbstractLinearCombo* self, char* (*morphismToString)(Morphism*)) {
    if (!self || !morphismToString) {
        return NULL;
    }

    if (self->term_count == 0) {
        char* zero_str = malloc(2);
        if (zero_str) strcpy(zero_str, "0");
        return zero_str;
    }

    StringBuilder sb;
    sb_init(&sb);
    for (int i = 0; i < self->term_count; i++) {
        Morphism* term = self->elements[i].term;
        BivariatePoly* coeff = self->elements[i].coeff;
        
        //generate strings
        char* coeff_str = bp_to_string(coeff); 
        char* term_str = morphismToString(term);

        //separator
        if (i > 0) {
            sb_append(&sb, " + ");
        }
        sb_append(&sb, "(");
        sb_append(&sb, coeff_str);
        sb_append(&sb, ") * ");
        sb_append(&sb, term_str);
        free(coeff_str); //freeing
        free(term_str);
    }
    return sb.str;
}

AbstractLinearCombo* AbstractLinearCombo_create(int initial_capacity) {
    AbstractLinearCombo* combo = (AbstractLinearCombo*)malloc(sizeof(AbstractLinearCombo));
    if (!combo) return NULL;
    combo->term_count = 0;
    combo->capacity = (initial_capacity > 0) ? initial_capacity : 4;
    combo->elements = (TermPair*)malloc(combo->capacity * sizeof(TermPair));
    if (!combo->elements) {
        free(combo);
        return NULL;
    }
    return combo;
}

void AbstractLinearCombo_addTerm(AbstractLinearCombo* self, Morphism* term, BivariatePoly* coeff) {
    if (!self) return;
    if (self->term_count >= self->capacity) {
        self->capacity = (self->capacity == 0) ? 4 : self->capacity * 2;
        TermPair* new_elements = (TermPair*)realloc(self->elements, self->capacity * sizeof(TermPair));
        if (!new_elements) {
            return;
        }
        self->elements = new_elements;
    }
    
    self->elements[self->term_count].term = term;
    self->elements[self->term_count].coeff = coeff;
    self->term_count++;
}

void AbstractLinearCombo_free(AbstractLinearCombo* self) {
    if (!self) return;
    if (self->elements) {
        for (int i = 0; i < self->term_count; i++) {
            if (self->elements[i].coeff) {
                bp_free(self->elements[i].coeff);
            }
        }
        free(self->elements);
    }
    free(self);
}
