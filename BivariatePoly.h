#ifndef BIVARIATE_POLY_H
#define BIVARIATE_POLY_H

typedef struct {
    int q_exp; // quantum degree
    int t_exp;// homological degree
    int coeff;
} Monomial;

typedef struct {
    Monomial* terms;// array of monomials
    int num_terms;// current number of terms
    int capacity;// allocated size
} BivariatePoly;
// create zero polynomial
BivariatePoly* bp_create(void);
// free memory
void bp_free(BivariatePoly* p);
// add 
BivariatePoly* bp_add(const BivariatePoly* p1, const BivariatePoly* p2);
// multiply
BivariatePoly* bp_multiply(const BivariatePoly* p1, const BivariatePoly* p2);
void bp_print(const BivariatePoly* p);

#endif