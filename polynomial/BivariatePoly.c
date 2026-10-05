#include "../KhMemory.h"
#include "BivariatePoly.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* ---- internal helpers ---- */

static int monomial_cmp(const Monomial *a, const Monomial *b)
{
    if (a->t_exp != b->t_exp)
        return (a->t_exp < b->t_exp) ? -1 : 1;
    if (a->q_exp != b->q_exp)
        return (a->q_exp < b->q_exp) ? -1 : 1;
    return 0;
}

static int monomial_cmp_qsort(const void *x, const void *y)
{
    return monomial_cmp((const Monomial *)x, (const Monomial *)y);
}

/* Ensure capacity for at least one more term; allocation failure is fatal. */
static void bp_ensure_capacity(BivariatePoly *p)
{
    if (p->num_terms < p->capacity)
        return;
    int new_cap = (p->capacity == 0) ? 8 : kh_int((int64_t)p->capacity * 2);
    Monomial *nt = (Monomial *)kh_realloc(p->terms, (size_t)new_cap * sizeof(Monomial));
    p->terms = nt;
    p->capacity = new_cap;
}

/* Integer power t0^e. t0 is a small integer (typically ±1); e may be negative
 * only when |t0|==1 (then result is still ±1). For |t0|!=1 and e<0 returns 0. */
static int int_pow(int base, int exp)
{
    if (exp == 0)
        return 1;
    if (base == 0)
        return 0;
    if (base == 1)
        return 1;
    if (base == -1)
        return (exp % 2 == 0) ? 1 : -1;
    if (exp < 0)
        return 0; /* cannot represent fractional values in Z coeffs */

    int result = 1;
    int b = base;
    int e = exp;
    while (e > 0) {
        if (e & 1)
            result = kh_int((int64_t)result * b);
        e >>= 1;
        if (e) b = kh_int((int64_t)b * b);
    }
    return result;
}

/* ---- lifecycle ---- */

BivariatePoly *bp_create(void)
{
    BivariatePoly *p = (BivariatePoly *)kh_malloc(sizeof(BivariatePoly));
    p->capacity = 8;
    p->num_terms = 0;
    p->terms = (Monomial *)kh_malloc((size_t)p->capacity * sizeof(Monomial));

    return p;
}

void bp_free(BivariatePoly *p)
{
    if (p) {
        free(p->terms);
        free(p);
    }
}

BivariatePoly *bp_copy(const BivariatePoly *p)
{
    if (!p)
        return bp_create();
    BivariatePoly *out = bp_create();

    for (int i = 0; i < p->num_terms; i++)
        bp_add_term(out, p->terms[i].q_exp, p->terms[i].t_exp, p->terms[i].coeff);
    return out;
}

/* ---- status ---- */

bool bp_is_zero(const BivariatePoly *p)
{
    return p == NULL || p->num_terms == 0;
}

bool bp_equals(const BivariatePoly *p1, const BivariatePoly *p2)
{
    if (p1 == p2)
        return true;
    if (bp_is_zero(p1) && bp_is_zero(p2))
        return true;
    if (bp_is_zero(p1) || bp_is_zero(p2))
        return false;
    if (p1->num_terms != p2->num_terms)
        return false;

    /* Both sides are kept sorted by bp_add_term / bp_sort. */
    for (int i = 0; i < p1->num_terms; i++) {
        if (p1->terms[i].q_exp != p2->terms[i].q_exp ||
            p1->terms[i].t_exp != p2->terms[i].t_exp ||
            p1->terms[i].coeff != p2->terms[i].coeff)
            return false;
    }
    return true;
}

/* ---- mutation ---- */

void bp_sort(BivariatePoly *p)
{
    if (!p || p->num_terms <= 1)
        return;
    qsort(p->terms, (size_t)p->num_terms, sizeof(Monomial), monomial_cmp_qsort);
}

void bp_add_term(BivariatePoly *p, int q, int t, int c)
{
    if (!p || c == 0)
        return;

    /* Merge into an existing like term if present. */
    for (int i = 0; i < p->num_terms; i++) {
        if (p->terms[i].q_exp == q && p->terms[i].t_exp == t) {
            p->terms[i].coeff = kh_int((int64_t)p->terms[i].coeff + c);
            if (p->terms[i].coeff == 0) {
                memmove(p->terms + i, p->terms + i + 1,
                        (size_t)(p->num_terms - i - 1) * sizeof(Monomial));
                p->num_terms--;
            }
            return;
        }
    }

    /* Append new term then restore sorted order. */
    bp_ensure_capacity(p);

    p->terms[p->num_terms].q_exp = q;
    p->terms[p->num_terms].t_exp = t;
    p->terms[p->num_terms].coeff = c;
    p->num_terms++;
    bp_sort(p);
}

void bp_shift(BivariatePoly *p, int dq, int dt)
{
    if (!p || (dq == 0 && dt == 0))
        return;
    for (int i = 0; i < p->num_terms; i++) {
        p->terms[i].q_exp = kh_int((int64_t)p->terms[i].q_exp + dq);
        p->terms[i].t_exp = kh_int((int64_t)p->terms[i].t_exp + dt);
    }
    /* Relative order of t then q is preserved under uniform shift. */
}

/* ---- algebra ---- */

BivariatePoly *bp_add(const BivariatePoly *p1, const BivariatePoly *p2)
{
    BivariatePoly *result = bp_create();

    if (p1) {
        for (int i = 0; i < p1->num_terms; i++)
            bp_add_term(result, p1->terms[i].q_exp, p1->terms[i].t_exp, p1->terms[i].coeff);
    }
    if (p2) {
        for (int i = 0; i < p2->num_terms; i++)
            bp_add_term(result, p2->terms[i].q_exp, p2->terms[i].t_exp, p2->terms[i].coeff);
    }
    return result;
}

BivariatePoly *bp_multiply(const BivariatePoly *p1, const BivariatePoly *p2)
{
    BivariatePoly *result = bp_create();

    if (!p1 || !p2)
        return result;
    for (int i = 0; i < p1->num_terms; i++) {
        for (int j = 0; j < p2->num_terms; j++) {
            int q = kh_int((int64_t)p1->terms[i].q_exp + p2->terms[j].q_exp);
            int t = kh_int((int64_t)p1->terms[i].t_exp + p2->terms[j].t_exp);
            int c = kh_int((int64_t)p1->terms[i].coeff * p2->terms[j].coeff);
            bp_add_term(result, q, t, c);
        }
    }
    return result;
}

BivariatePoly *bp_eval_t(const BivariatePoly *p, int t0)
{
    BivariatePoly *result = bp_create();
    if (!p)
        return result;

    for (int i = 0; i < p->num_terms; i++) {
        int weight = int_pow(t0, p->terms[i].t_exp);
        int c = kh_int((int64_t)p->terms[i].coeff * weight);
        /* Collapse onto the q-axis (t_exp = 0). */
        bp_add_term(result, p->terms[i].q_exp, 0, c);
    }
    return result;
}

BivariatePoly *bp_from_ranks(int **betti,
                             int h_min, int h_max,
                             int q_min, int q_max)
{
    BivariatePoly *result = bp_create();
    if (!betti)
        return result;
    if (h_max < h_min || q_max < q_min)
        return result;

    int q_count = q_max - q_min + 1;
    for (int h = h_min; h <= h_max; h++) {
        int hi = h - h_min;
        if (!betti[hi])
            continue;
        for (int q = q_min; q <= q_max; q++) {
            int qi = q - q_min;
            if (qi < 0 || qi >= q_count)
                continue;
            int rank = betti[hi][qi];
            if (rank > 0)
                bp_add_term(result, q, h, rank);
        }
    }
    return result;
}

/* ---- output ---- */

void bp_print(const BivariatePoly *p)
{
    if (!p || p->num_terms == 0) {
        printf("0\n");
        return;
    }

    for (int i = 0; i < p->num_terms; i++) {
        int coeff = p->terms[i].coeff;
        if (coeff == 0)
            continue;

        if (i == 0) {
            if (coeff < 0)
                printf("-");
        } else {
            if (coeff < 0)
                printf(" - ");
            else
                printf(" + ");
        }

        int64_t c = coeff < 0 ? -(int64_t)coeff : coeff;
        int q = p->terms[i].q_exp;
        int t = p->terms[i].t_exp;
        int has_vars = (q != 0) || (t != 0);

        if (c != 1 || !has_vars)
            printf("%lld", (long long)c);

        if (q != 0) {
            if (q == 1)
                printf("q");
            else if (q == -1)
                printf("q^{-1}");
            else
                printf("q^{%d}", q);
        }
        if (t != 0) {
            if (t == 1)
                printf("t");
            else if (t == -1)
                printf("t^{-1}");
            else
                printf("t^{%d}", t);
        }
    }
    printf("\n");
}

char *bp_to_string(const BivariatePoly *p)
{
    if (!p || p->num_terms == 0) {
        char *zero_str = (char *)kh_malloc(2);

        zero_str[0] = '0';
        zero_str[1] = '\0';
        return zero_str;
    }

    size_t buf_size = kh_size_mul((size_t)p->num_terms + 4, 64);
    char *buffer = (char *)kh_malloc(buf_size);

    size_t pos = 0;
    buffer[0] = '\0';

    for (int i = 0; i < p->num_terms; i++) {
        int coeff = p->terms[i].coeff;
        if (coeff == 0)
            continue;

        if (pos == 0) {
            if (coeff < 0) {
                buffer[pos++] = '-';
                buffer[pos] = '\0';
            }
        } else {
            int written;
            if (coeff < 0)
                written = snprintf(buffer + pos, buf_size - pos, " - ");
            else
                written = snprintf(buffer + pos, buf_size - pos, " + ");
            if (written < 0 || (size_t)written >= buf_size - pos) {
                buffer[buf_size - 1] = '\0';
                return buffer;
            }
            pos += (size_t)written;
        }

        int64_t c = coeff < 0 ? -(int64_t)coeff : coeff;
        int q = p->terms[i].q_exp;
        int t = p->terms[i].t_exp;
        int has_vars = (q != 0) || (t != 0);
        int written;

        if (c != 1 || !has_vars) {
            written = snprintf(buffer + pos, buf_size - pos, "%lld", (long long)c);
            if (written < 0 || (size_t)written >= buf_size - pos) {
                buffer[buf_size - 1] = '\0';
                return buffer;
            }
            pos += (size_t)written;
        }
        if (q != 0) {
            if (q == 1)
                written = snprintf(buffer + pos, buf_size - pos, "q");
            else if (q == -1)
                written = snprintf(buffer + pos, buf_size - pos, "q^{-1}");
            else
                written = snprintf(buffer + pos, buf_size - pos, "q^{%d}", q);
            if (written < 0 || (size_t)written >= buf_size - pos) {
                buffer[buf_size - 1] = '\0';
                return buffer;
            }
            pos += (size_t)written;
        }
        if (t != 0) {
            if (t == 1)
                written = snprintf(buffer + pos, buf_size - pos, "t");
            else if (t == -1)
                written = snprintf(buffer + pos, buf_size - pos, "t^{-1}");
            else
                written = snprintf(buffer + pos, buf_size - pos, "t^{%d}", t);
            if (written < 0 || (size_t)written >= buf_size - pos) {
                buffer[buf_size - 1] = '\0';
                return buffer;
            }
            pos += (size_t)written;
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
