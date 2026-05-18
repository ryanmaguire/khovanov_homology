#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../polynomial/polynomial/BivariatePoly.h"
#include "../Morphism.h"
#include "LinearComboMap.h"

typedef struct Morphism {
    const char* name;
} TestMorphism;

static char* mock_morphismToString(Morphism* m) {
    if (!m) return NULL;
    TestMorphism* tm = (TestMorphism*)m;
    if (!tm->name) return NULL;
    char* str = malloc(strlen(tm->name) + 1);
    if (!str) return NULL;
    strcpy(str, tm->name);
    return str;
}

BivariatePoly* create_poly(int c, int q, int t) {
    BivariatePoly* p = bp_create();
    if (c != 0) {
        p->terms[0].coeff = c;
        p->terms[0].q_exp = q;
        p->terms[0].t_exp = t;
        p->num_terms = 1;
    }
    return p;
}

int main(void) {
    TestMorphism m1 = { .name = "Morphism_A" };
    TestMorphism m2 = { .name = "Morphism_B" };

    LinearComboMap* myCombo = LinearComboMap_create(4);
    if (!myCombo) return 1;

    BivariatePoly* p1 = create_poly(2, 1, 0); 
    LinearComboMap_putOrAdd(myCombo, (Morphism*)&m1, p1);
    
    char* str1 = AbstractLinearCombo_toString(myCombo->combo, mock_morphismToString);
    printf("Add: %s\n", str1 ? str1 : "NULL");
    free(str1);

    BivariatePoly* p2 = create_poly(3, 1, 0);
    LinearComboMap_putOrAdd(myCombo, (Morphism*)&m1, p2);

    char* str2 = AbstractLinearCombo_toString(myCombo->combo, mock_morphismToString);
    printf("Merge: %s\n", str2 ? str2 : "NULL");
    free(str2);

    BivariatePoly* p3 = create_poly(-1, 0, 2);
    LinearComboMap_putOrAdd(myCombo, (Morphism*)&m2, p3);

    char* str3 = AbstractLinearCombo_toString(myCombo->combo, mock_morphismToString);
    printf("New: %s\n", str3 ? str3 : "NULL");
    free(str3);

    BivariatePoly* scalar = create_poly(2, 1, 0);
    LinearComboMap_multiply(myCombo, scalar); 

    char* str4 = AbstractLinearCombo_toString(myCombo->combo, mock_morphismToString);
    printf("Multiply: %s\n", str4 ? str4 : "NULL");
    free(str4);
    bp_free(scalar); 

    AbstractLinearCombo* emptyCombo = AbstractLinearCombo_create(0);
    int zeroTerms = AbstractLinearCombo_numberOfTerms(emptyCombo);
    bool isZero = AbstractLinearCombo_isZero(emptyCombo);
    
    if (zeroTerms == 0 && isZero == true) {
        printf("Zero: OK\n");
    }
    AbstractLinearCombo_free(emptyCombo);

    LinearComboMap_free(myCombo);

    printf("Done.\n");
    return 0;
}