#include "BivariatePoly.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// helper: create empty poly with space for 8 terms
BivariatePoly* bp_create(void) {
    BivariatePoly* p = malloc(sizeof(BivariatePoly));
    p->capacity = 8;
    p->num_terms = 0;
    p->terms = malloc(p->capacity * sizeof(Monomial));
    return p;
}
void bp_free(BivariatePoly* p) {  
    if (p) {
        free(p->terms);
        free(p);
    }
}
// add term, merge like terms
static void bp_add_term(BivariatePoly* p, int q, int t, int c) {
    if (c == 0) return;
    for (int i = 0; i < p->num_terms; i++) {
        if (p->terms[i].q_exp == q && p->terms[i].t_exp == t) {
            p->terms[i].coeff += c;
            if (p->terms[i].coeff == 0) {  // remove zero
                memmove(p->terms + i, p->terms + i + 1,
                        (p->num_terms - i - 1) * sizeof(Monomial));
                p->num_terms--;
            }
            return;
        }
    }
    // new term
    if (p->num_terms == p->capacity) {
        p->capacity *= 2;
        p->terms = realloc(p->terms, p->capacity * sizeof(Monomial));
    }
    p->terms[p->num_terms].q_exp = q;
    p->terms[p->num_terms].t_exp = t;
    p->terms[p->num_terms].coeff = c;
    p->num_terms++;
}
BivariatePoly* bp_add(const BivariatePoly* p1, const BivariatePoly* p2) {
    BivariatePoly* result = bp_create();
    for (int i = 0; i < p1->num_terms; i++)
        bp_add_term(result, p1->terms[i].q_exp, p1->terms[i].t_exp, p1->terms[i].coeff);
    for (int i = 0; i < p2->num_terms; i++)
        bp_add_term(result, p2->terms[i].q_exp, p2->terms[i].t_exp, p2->terms[i].coeff);
    return result;
}
BivariatePoly* bp_multiply(const BivariatePoly* p1, const BivariatePoly* p2) {
    BivariatePoly* result = bp_create();
    for (int i = 0; i < p1->num_terms; i++) {
        for (int j = 0; j < p2->num_terms; j++) {
            int q = p1->terms[i].q_exp + p2->terms[j].q_exp;
            int t = p1->terms[i].t_exp + p2->terms[j].t_exp;
            int c = p1->terms[i].coeff * p2->terms[j].coeff;
            bp_add_term(result, q, t, c);
        }
    }
    return result;
}
void bp_print(const BivariatePoly* p) {  //output
    if (p->num_terms == 0) {
        printf("0");
        return;
    }
    for (int i = 0; i < p->num_terms; i++) {
        if (i > 0 && p->terms[i].coeff > 0) printf(" + ");
        else if (p->terms[i].coeff < 0) printf(" - ");
        int c = abs(p->terms[i].coeff);
        if (c != 1 || (p->terms[i].q_exp == 0 && p->terms[i].t_exp == 0))
            printf("%d", c);
        if (p->terms[i].q_exp != 0) printf("q^%d", p->terms[i].q_exp);
        if (p->terms[i].t_exp != 0) printf("t^%d", p->terms[i].t_exp);
    }
    printf("\n");
}

char* bp_to_string(const BivariatePoly* p) {
    if (!p || p->num_terms == 0) {
        char* zero_str = malloc(2);
        if (!zero_str) return NULL;
        zero_str[0] = '0';
        zero_str[1] = '\0';
        return zero_str;
    }

    size_t buf_size = 256 + (size_t)p->num_terms * 64;
    char* buffer = malloc(buf_size);
    if (!buffer) return NULL;
    size_t pos = 0;

    for (int i = 0; i < p->num_terms; i++) {
        int coeff = p->terms[i].coeff;
        if (coeff == 0) continue;

        if (pos == 0) {
            if (coeff < 0) {
                buffer[pos++] = '-';
                buffer[pos] = '\0';
            }
        } else {
            if (coeff < 0) {
                pos += (size_t)snprintf(buffer + pos, buf_size - pos, " - ");
            } else {
                pos += (size_t)snprintf(buffer + pos, buf_size - pos, " + ");
            }
        }

        int c = abs(coeff);
        if (c != 1 || (p->terms[i].q_exp == 0 && p->terms[i].t_exp == 0)) {
            pos += (size_t)snprintf(buffer + pos, buf_size - pos, "%d", c);
        }
        if (p->terms[i].q_exp != 0) {
            pos += (size_t)snprintf(buffer + pos, buf_size - pos, "q^%d", p->terms[i].q_exp);
        }
        if (p->terms[i].t_exp != 0) {
            pos += (size_t)snprintf(buffer + pos, buf_size - pos, "t^%d", p->terms[i].t_exp);
        }

        if (pos >= buf_size) {
            buffer[buf_size - 1] = '\0';
            break;
        }
    }

    if (pos == 0) {
        buffer[0] = '0';
        buffer[1] = '\0';
    }

    return buffer;
}
