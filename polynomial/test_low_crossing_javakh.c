/*
 * Unit checks for 4-, 5-, 6-crossing knots (Ryan task 3).
 *
 * What this does
 * --------------
 * 1. Build standard Rolfsen knots from Knot Atlas PD codes (1-based → 0-based).
 * 2. Run kh_poincare() (true free Poincaré path).
 * 3. Mandatory: Kh(q,-1) == (q + q^{-1}) * V(q^2) using Knot Atlas Jones.
 * 4. For 5_1: compare free Poincaré to Knot Atlas / JavaKh table (hard gold).
 * 5. Print full Kh for every knot so you can diff against JavaKh by hand.
 *
 * Knots covered: 4_1, 5_1, 5_2, 6_1 (one per crossing number 4–6, plus extra 5_2).
 *
 * Build (from repo root, after steps 1–8 files are in place):
 *   gcc -O2 -o test_low_crossing \
 *     polynomial/test/test_low_crossing_javakh.c \
 *     polynomial/homology/KhPoincare.c \
 *     polynomial/polynomial/BivariatePoly.c \
 *     IntegerMatrix.c \
 *     -I. -Ipolynomial/homology -Ipolynomial/polynomial -Ipolynomial/encoder
 *
 *   ./test_low_crossing
 *
 * How to finish JavaKh parity
 * ---------------------------
 * Run JavaKh on the same PD, copy its Kh[q,t] into the gold_* helpers below,
 * then switch the corresponding expect_poly from "print only" to hard assert.
 */

#include "KhPoincare.h"
#include "BivariatePoly.h"
#include "KnotEncoder.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_failures = 0;

static void expect_poly(const char *name, const BivariatePoly *got, const BivariatePoly *want)
{
    if (bp_equals(got, want)) {
        printf("  PASS  %s\n", name);
        printf("        ");
        bp_print(got);
    } else {
        printf("  FAIL  %s\n", name);
        printf("        got:  ");
        bp_print(got);
        printf("        want: ");
        bp_print(want);
        g_failures++;
    }
}

/* ---- PD helpers: Knot Atlas X[a,b,c,d] is 1-based; knot_from_pd converts ---- */

static Knot *knot_from_pd_rows(int rows[][4], int m)
{
    return knot_from_pd(rows, m);
}

/* 4_1 figure-eight: PD X[4,2,5,1] X[8,6,1,5] X[6,3,7,4] X[2,7,3,8] */
static Knot *knot_4_1(void)
{
    int pd[4][4] = {
        {4, 2, 5, 1},
        {8, 6, 1, 5},
        {6, 3, 7, 4},
        {2, 7, 3, 8}
    };
    return knot_from_pd_rows(pd, 4);
}

/* 5_1 cinquefoil / T(2,5): X[1,6,2,7] X[3,8,4,9] X[5,10,6,1] X[7,2,8,3] X[9,4,10,5] */
static Knot *knot_5_1(void)
{
    int pd[5][4] = {
        {1, 6, 2, 7},
        {3, 8, 4, 9},
        {5, 10, 6, 1},
        {7, 2, 8, 3},
        {9, 4, 10, 5}
    };
    return knot_from_pd_rows(pd, 5);
}

/* 5_2 three-twist: X[1,4,2,5] X[3,8,4,9] X[5,10,6,1] X[9,6,10,7] X[7,2,8,3] */
static Knot *knot_5_2(void)
{
    int pd[5][4] = {
        {1, 4, 2, 5},
        {3, 8, 4, 9},
        {5, 10, 6, 1},
        {9, 6, 10, 7},
        {7, 2, 8, 3}
    };
    return knot_from_pd_rows(pd, 5);
}

/* 6_1 stevedore: Knot Atlas X[1,4,2,5] X[3,8,4,9] X[5,10,6,11] X[9,12,10,7] X[11,2,12,3] X[7,6,8,1]
 * (standard Rolfsen 6_1 PD; verify writhe/Jones if chi fails) */
static Knot *knot_6_1(void)
{
    int pd[6][4] = {
        {1, 4, 2, 5},
        {3, 8, 4, 9},
        {5, 10, 6, 11},
        {9, 12, 10, 7},
        {11, 2, 12, 3},
        {7, 6, 8, 1}
    };
    return knot_from_pd_rows(pd, 6);
}

/* ---- Jones χ_q = (q + q^{-1}) V(q^2) from Knot Atlas Jones ---- */

/* 4_1: V = q^2 + q^{-2} - q - q^{-1} + 1
 * V(q^2) = q^4 + q^{-4} - q^2 - q^{-2} + 1
 * χ = (q+q^{-1}) V(q^2) = q^5 + q - q^{-1} + q^{-5}
 */
static BivariatePoly *jones_chi_4_1(void)
{
    BivariatePoly *p = bp_create();
    bp_add_term(p, 5, 0, 1);
    bp_add_term(p, 1, 0, 1);
    bp_add_term(p, -1, 0, -1);
    bp_add_term(p, -5, 0, 1);
    return p;
}

/* 5_1: V = -q^{-7} + q^{-6} - q^{-5} + q^{-4} + q^{-2}
 * (Knot Atlas). χ = (q+q^{-1}) V(q^2).
 * Atlas also lists Kh(q,-1) = -q^{-15} + q^{-7} + q^{-5} + q^{-3}
 * (from the Khovanov Homology page verification).
 */
static BivariatePoly *jones_chi_5_1(void)
{
    BivariatePoly *p = bp_create();
    bp_add_term(p, -15, 0, -1);
    bp_add_term(p, -7, 0, 1);
    bp_add_term(p, -5, 0, 1);
    bp_add_term(p, -3, 0, 1);
    return p;
}

/* 5_2: V = -q^{-6} + q^{-5} - q^{-4} + 2 q^{-3} - q^{-2} + q^{-1}
 * χ expanded from (q+q^{-1}) V(q^2):
 * V(q^2) = -q^{-12} + q^{-10} - q^{-8} + 2 q^{-6} - q^{-4} + q^{-2}
 * χ = -q^{-13} - q^{-11} + q^{-9} + q^{-7} + q^{-5} - q^{-3} + q^{-1}
 *     + (-q^{-11} + q^{-9} - q^{-7} + 2 q^{-5} - q^{-3} + q^{-1}) wait — compute carefully:
 * (q + q^{-1}) * V(q^2) =
 * q*(-q^{-12}+q^{-10}-q^{-8}+2q^{-6}-q^{-4}+q^{-2})
 * + q^{-1}*(-q^{-12}+q^{-10}-q^{-8}+2q^{-6}-q^{-4}+q^{-2})
 * = -q^{-11}+q^{-9}-q^{-7}+2q^{-5}-q^{-3}+q^{-1}
 *   -q^{-13}+q^{-11}-q^{-9}+2q^{-7}-q^{-5}+q^{-3}
 * = -q^{-13} + ( -1+1 ) q^{-11} + (-1+1) q^{-9} + (-1+2) q^{-7} + (2-1) q^{-5} + (-1+1) q^{-3} + q^{-1}
 * = -q^{-13} + q^{-7} + q^{-5} + q^{-1}
 */
static BivariatePoly *jones_chi_5_2(void)
{
    BivariatePoly *p = bp_create();
    bp_add_term(p, -13, 0, -1);
    bp_add_term(p, -7, 0, 1);
    bp_add_term(p, -5, 0, 1);
    bp_add_term(p, -1, 0, 1);
    return p;
}

/* 6_1: V from Knot Atlas is typically
 * V = q^{-2} - q^{-1} + 1 - q + q^2 - q^3 + q^4  (one orientation)
 * We compute chi from the engine and only soft-check by printing unless
 * you fill the gold after JavaKh. For a hard Jones check, fill jones_chi_6_1
 * from Atlas then enable the assert below.
 */
static BivariatePoly *jones_chi_6_1_placeholder(void)
{
    /* Empty marker: Jones hard-check skipped until filled from Atlas/JavaKh. */
    return bp_create();
}

/* Knot Atlas free Poincaré for 5_1 (characteristic 0):
 * q^{-5} + q^{-3} + q^{-15} t^{-5} + q^{-11} t^{-4} + q^{-11} t^{-3} + q^{-7} t^{-2}
 */
static BivariatePoly *gold_5_1(void)
{
    BivariatePoly *p = bp_create();
    bp_add_term(p, -5, 0, 1);
    bp_add_term(p, -3, 0, 1);
    bp_add_term(p, -15, -5, 1);
    bp_add_term(p, -11, -4, 1);
    bp_add_term(p, -11, -3, 1);
    bp_add_term(p, -7, -2, 1);
    return p;
}

static void run_one(const char *name, Knot *k, BivariatePoly *want_chi, BivariatePoly *want_kh)
{
    printf("[%s]  m=%d writhe=%d\n", name, k->m, k->writhe);
    BivariatePoly *kh = kh_poincare(k);
    printf("  Kh(q,t) = ");
    bp_print(kh);

    if (want_chi && want_chi->num_terms > 0) {
        BivariatePoly *chi = bp_eval_t(kh, -1);
        char buf[128];
        snprintf(buf, sizeof(buf), "%s  Kh(q,-1) == (q+q^{-1})V(q^2)", name);
        expect_poly(buf, chi, want_chi);
        bp_free(chi);
    } else {
        printf("  SKIP  Jones chi hard-check (fill gold from JavaKh/Atlas)\n");
        BivariatePoly *chi = bp_eval_t(kh, -1);
        printf("  Kh(q,-1) = ");
        bp_print(chi);
        bp_free(chi);
    }

    if (want_kh) {
        char buf[128];
        snprintf(buf, sizeof(buf), "%s  free Poincaré vs JavaKh/Atlas gold", name);
        expect_poly(buf, kh, want_kh);
    } else {
        printf("  NOTE  paste this Kh into JavaKh diff; then hardcode gold_%s\n", name);
    }

    bp_free(kh);
    if (want_chi)
        bp_free(want_chi);
    if (want_kh)
        bp_free(want_kh);
    knot_free(k);
    printf("\n");
}

int main(void)
{
    printf("=== Low-crossing unit checks (4/5/6) vs Jones + JavaKh gold ===\n\n");

    run_one("4_1", knot_4_1(), jones_chi_4_1(), NULL);
    run_one("5_1", knot_5_1(), jones_chi_5_1(), gold_5_1());
    run_one("5_2", knot_5_2(), jones_chi_5_2(), NULL);
    run_one("6_1", knot_6_1(), jones_chi_6_1_placeholder(), NULL);

    if (g_failures == 0) {
        printf("=== ALL HARD CHECKS PASSED ===\n");
        printf("Next: run JavaKh on the same PDs and fill remaining gold_* polynomials.\n");
        return 0;
    }
    printf("=== %d FAILURE(S) ===\n", g_failures);
    return 1;
}
