#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LinearComboMap.h"
#include "ZeroLinearCombo.h"
#include "../../polynomial/polynomial/BivariatePoly.h"
#include "../Morphism.h"

typedef struct {
    Morphism base;
    const char* name;
} TestMorphism;

static char* mock_morphismToString(Morphism* m) {
    if (!m) return NULL;
    const TestMorphism* tm = (const TestMorphism*)m;
    if (!tm->name) return NULL;
    char* str = malloc(strlen(tm->name) + 1);
    if (!str) return NULL;
    strcpy(str, tm->name);
    return str;
}

static int mock_morphismCompareTo(Morphism* self, Morphism* other) {
    if (self == other) return 0;
    if (!self) return -1;
    if (!other) return 1;
    const TestMorphism* a = (const TestMorphism*)self;
    const TestMorphism* b = (const TestMorphism*)other;
    if (a->name == b->name) return 0;
    if (!a->name) return -1;
    if (!b->name) return 1;
    return strcmp(a->name, b->name);
}

//create a single-term polynomial for our tests
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
    printf("--- Khovanov Homology Architecture Test ---\n\n");

    TestMorphism m1 = {0};
    m1.base.compareTo = mock_morphismCompareTo;
    m1.name = "Morphism_A";
    TestMorphism m2 = {0};
    m2.base.compareTo = mock_morphismCompareTo;
    m2.name = "Morphism_B";

    LinearComboMap myCombo;
    LinearComboMap_init(&myCombo);
    
    myCombo.base.morphismToString = mock_morphismToString; 
    printf("LinearComboMap initialized.\n");

    BivariatePoly* p1 = create_poly(2, 1, 0); 
    LinearComboMap_addTerm(&myCombo, (Morphism*)&m1, p1);
    
    char* str1 = AbstractLinearCombo_toString((AbstractLinearCombo*)&myCombo);
    printf("Test 1 (Add): %s\n", str1 ? str1 : "NULL");
    free(str1);

    // erge a term
    //add (3q^1 * m1) should merge to become (5q^1 * m1)
    BivariatePoly* p2 = create_poly(3, 1, 0);
    LinearComboMap_addTerm(&myCombo, (Morphism*)&m1, p2);

    char* str2 = AbstractLinearCombo_toString((AbstractLinearCombo*)&myCombo);
    printf("Test 2 (Merge): %s\n", str2 ? str2 : "NULL");
    free(str2);

    BivariatePoly* p3 = create_poly(-1, 0, 2);
    LinearComboMap_addTerm(&myCombo, (Morphism*)&m2, p3);

    char* str3 = AbstractLinearCombo_toString((AbstractLinearCombo*)&myCombo);
    printf("Test 3 (Add New): %s\n", str3 ? str3 : "NULL");
    free(str3);

    //multiply the entire combo by (2q^1)
    BivariatePoly* scalar = create_poly(2, 1, 0);
    LinearComboMap_multiply(&myCombo, scalar); 

    char* str4 = AbstractLinearCombo_toString((AbstractLinearCombo*)&myCombo);
    printf("Test 4 (Multiply): %s\n", str4 ? str4 : "NULL");
    free(str4);
    bp_free(scalar); 

    //ZeroLinearCombo
    int zeroTerms = ZeroLinearCombo_numberOfTerms(NULL);
    bool isZero = ZeroLinearCombo_isZero(NULL);
    
    if (zeroTerms == 0 && isZero == true) {
        printf("\n[SUCCESS] ZeroLinearCombo reported 0 terms and isZero = true.\n");
    }
    for (int i = 0; i < myCombo.size; i++) {
        bp_free(myCombo.terms[i].coefficient);
    }
    free(myCombo.terms);

    printf("All tests complete.\n");

    return 0;
}
