#include "BivariatePoly.h"
#include <stdio.h>

int main() {
    // Test 1
    BivariatePoly* p1 = bp_create();
    bp_add_term(p1, 1, 0, 1);   //q
    bp_add_term(p1, 0, 1, 1);   //t
    BivariatePoly* p2 = bp_create();
    bp_add_term(p2, 1, 0, 1);   //q
    bp_add_term(p2, 0, 0, 2);   //2
    BivariatePoly* sum = bp_add(p1, p2);
    printf("p1+p2= ");
    bp_print(sum);
    // Test 2:
    BivariatePoly* prod = bp_multiply(p1, p2);
    printf("p1 * p2 = ");
    bp_print(prod);
    // Test 3: trefoil example

    BivariatePoly* trefoil = bp_create();
    bp_add_term(trefoil, -2, 0, 1);
    bp_add_term(trefoil, 0, -1, 1);
    bp_add_term(trefoil, 2, -2, 1);
    bp_add_term(trefoil, 4, -1, 1);
    bp_add_term(trefoil, 6, 0, 1);
    printf("Trefoil Khovanov polynomial (manual check): ");
    bp_print(trefoil);

    bp_free(p1); bp_free(p2); bp_free(sum); bp_free(prod); bp_free(trefoil);
    return 0;
}
