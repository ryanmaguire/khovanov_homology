#include "AbstractLinearCombo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char* str;
    size_t length;
    size_t capacity;
} StringBuilder;

static void sb_init(StringBuilder* sb) {
    sb->capacity = 64;
    sb->str = malloc(sb->capacity);
    if (sb->str) {
        sb->str[0] = '\0';
    }
    sb->length = 0;
}

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
Morphism* AbstractLinearCombo_firstTerm(AbstractLinearCombo* self) {
    if (!self || !self->terms) return NULL;
    
    MorphismCollection* terms = self->terms(self);
    if (!terms || !terms->iterator) return NULL;

    MorphismIterator* it = terms->iterator(terms);
    Morphism* first = NULL;
    
    if (it && it->hasNext(it)) {
        first = it->next(it);
    }
    
    if (it && it->free) it->free(it); 
    return first;
}

Ring* AbstractLinearCombo_firstCoefficient(AbstractLinearCombo* self) {
    Morphism* first = AbstractLinearCombo_firstTerm(self);
    if (first != NULL && self->getCoefficient) {
        return self->getCoefficient(self, first);
    }
    return NULL;
}

bool AbstractLinearCombo_isZero(AbstractLinearCombo* self) {
    if (!self || !self->terms) return true;
    MorphismCollection* terms = self->terms(self);
    return terms->isEmpty(terms);
}

int AbstractLinearCombo_numberOfTerms(AbstractLinearCombo* self) {
    if (!self || !self->terms) return 0;
    MorphismCollection* terms = self->terms(self);
    return terms->size(terms);
}

char* AbstractLinearCombo_toString(AbstractLinearCombo* self) {
    if (!self || !self->terms || !self->ringToString || !self->morphismToString) {
        return NULL;
    }

    StringBuilder sb;
    sb_init(&sb);

    MorphismCollection* terms_col = self->terms(self);
    MorphismIterator* i = terms_col->iterator(terms_col);

    if (i && i->hasNext(i)) {
        Morphism* term = i->next(i);
        Ring* coeff = self->getCoefficient(self, term);

        char* coeff_str = self->ringToString(coeff);
        char* term_str = self->morphismToString(term);

        sb_append(&sb, coeff_str);
        sb_append(&sb, " * ");
        sb_append(&sb, term_str);

        free(coeff_str);
        free(term_str);
    }

    while (i && i->hasNext(i)) {
        Morphism* term = i->next(i);
        Ring* coeff = self->getCoefficient(self, term);

        char* coeff_str = self->ringToString(coeff);
        char* term_str = self->morphismToString(term);

        sb_append(&sb, " + ");
        sb_append(&sb, coeff_str);
        sb_append(&sb, " * ");
        sb_append(&sb, term_str);

        free(coeff_str);
        free(term_str);
    }

    if (i && i->free) i->free(i);

    return sb.str;
}