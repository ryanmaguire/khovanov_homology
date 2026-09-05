/*
 * Step 5 — gold-table regression for Khovanov Poincaré polynomials.
 *
 * Convention (interim report §3.10–3.11, Knot Atlas):
 *   Kh(L)(q,t) = Σ t^r q^j dim Kh^{r,j}(L)
 *   Kh(q, -1)  = χ_q = (q + q^{-1}) V_L(q^2)
 *
 * Gold free ranks (rational Poincaré):
 *   Unknot:        q^{-1} + q
 *   Right trefoil: q + q^3 + q^5 t^2 + q^9 t^3
 *   Left trefoil:  q^{-1} + q^{-3} + q^{-5} t^{-2} + q^{-9} t^{-3}
 *     (mirror of right: (q,t) ↦ (q^{-1}, t^{-1}))
 *
 * Jones polynomials V(q) used for the self-check (Wikipedia / report):
 *   Unknot:        V = 1
 *   Right trefoil: V = q^{-1} + q^{-3} - q^{-4}
 *     so (q+q^{-1}) V(q^2) = q + q^3 + q^5 - q^9
 *   Left trefoil:  V = q + q^3 - q^4
 *     so (q+q^{-1}) V(q^2) = q^{-1} + q^{-3} + q^{-5} - q^{-9}
 *
 * Build (repo root, after steps 1–4 files are in place):
 *   gcc -O2 -o test_gold_khovanov \
 *     polynomial/test/test_gold_khovanov.c \
 *     polynomial/homology/KhPoincare.c \
 *     polynomial/polynomial/BivariatePoly.c \
 *     IntegerMatrix.c \
 *     -I. -Ipolynomial/homology -Ipolynomial/polynomial -Ipolynomial/encoder
 */

#include "../homology/KhPoincare.h"
#include "../polynomial/BivariatePoly.h"
#include "../encoder/KnotEncoder.h"

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

/* χ_q expected = (q + q^{-1}) * V(q^2) assembled as a univariate in q (t=0). */
static BivariatePoly *jones_chi_unknot(void)
{
    /* V=1 → (q + q^{-1}) * 1 */
    BivariatePoly *p = bp_create();
    bp_add_term(p, 1, 0, 1);
    bp_add_term(p, -1, 0, 1);
    return p;
}

static BivariatePoly *jones_chi_right_trefoil(void)
{
    /* (q+q^{-1})(q^2 + q^6 - q^8) = q + q^3 + q^5 - q^9 */
    BivariatePoly *p = bp_create();
    bp_add_term(p, 1, 0, 1);
    bp_add_term(p, 3, 0, 1);
    bp_add_term(p, 5, 0, 1);
    bp_add_term(p, 9, 0, -1);
    return p;
}

static BivariatePoly *jones_chi_left_trefoil(void)
{
    /* mirror: q^{-1} + q^{-3} + q^{-5} - q^{-9} */
    BivariatePoly *p = bp_create();
    bp_add_term(p, -1, 0, 1);
    bp_add_term(p, -3, 0, 1);
    bp_add_term(p, -5, 0, 1);
    bp_add_term(p, -9, 0, -1);
    return p;
}

static BivariatePoly *gold_unknot(void)
{
    BivariatePoly *p = bp_create();
    bp_add_term(p, -1, 0, 1);
    bp_add_term(p, 1, 0, 1);
    return p;
}

static BivariatePoly *gold_right_trefoil(void)
{
    BivariatePoly *p = bp_create();
    bp_add_term(p, 1, 0, 1);
    bp_add_term(p, 3, 0, 1);
    bp_add_term(p, 5, 2, 1);
    bp_add_term(p, 9, 3, 1);
    return p;
}

static BivariatePoly *gold_left_trefoil(void)
{
    BivariatePoly *p = bp_create();
    bp_add_term(p, -1, 0, 1);
    bp_add_term(p, -3, 0, 1);
    bp_add_term(p, -5, -2, 1);
    bp_add_term(p, -9, -3, 1);
    return p;
}

/* Right-handed trefoil: 0-based PD, writhe = +3
 * (verified to yield literature gold q + q^3 + q^5 t^2 + q^9 t^3) */
static Knot *knot_create_right_trefoil(void)
{
    static const int pd_data[3][4] = {
        {4, 0, 5, 1},
        {1, 3, 2, 4},
        {3, 5, 0, 2}
    };
    Knot *knot = (Knot *)malloc(sizeof(Knot));
    if (!knot)
        return NULL;
    knot->m = 3;
    knot->writhe = 3;
    knot->edges = NULL;
    knot->pd = (int **)malloc(3 * sizeof(int *));
    if (!knot->pd) {
        free(knot);
        return NULL;
    }
    for (int i = 0; i < 3; i++) {
        knot->pd[i] = (int *)malloc(4 * sizeof(int));
        if (!knot->pd[i]) {
            for (int j = 0; j < i; j++)
                free(knot->pd[j]);
            free(knot->pd);
            free(knot);
            return NULL;
        }
        for (int j = 0; j < 4; j++)
            knot->pd[i][j] = pd_data[i][j];
    }
    return knot;
}

/* Unknot: m = 0, no crossings — one unknotted circle */
static Knot *knot_create_unknot(void)
{
    Knot *knot = (Knot *)malloc(sizeof(Knot));
    if (!knot)
        return NULL;
    knot->m = 0;
    knot->writhe = 0;
    knot->edges = NULL;
    knot->pd = NULL;
    return knot;
}

static void check_jones(const char *name, const BivariatePoly *kh, BivariatePoly *want_chi)
{
    BivariatePoly *chi = bp_eval_t(kh, -1);
    char buf[128];
    snprintf(buf, sizeof(buf), "%s  Kh(q,-1) == (q+q^{-1})V(q^2)", name);
    expect_poly(buf, chi, want_chi);
    bp_free(chi);
    bp_free(want_chi);
}

int main(void)
{
    printf("=== Gold table: Khovanov Poincaré ===\n\n");

    /* ---- Unknot ---- */
    {
        printf("[unknot]\n");
        Knot *u = knot_create_unknot();
        BivariatePoly *kh = kh_poincare(u);
        BivariatePoly *want = gold_unknot();
        expect_poly("unknot Kh(q,t)", kh, want);
        check_jones("unknot", kh, jones_chi_unknot());
        bp_free(want);
        bp_free(kh);
        knot_free(u);
        printf("\n");
    }

    /* ---- Left trefoil ---- */
    {
        printf("[left trefoil]\n");
        Knot *k = knot_create_left_trefoil();
        KhTorsionReport tors;
        kh_torsion_report_init(&tors);
        BivariatePoly *kh = kh_poincare_with_torsion(k, &tors);
        BivariatePoly *want = gold_left_trefoil();
        expect_poly("left trefoil Kh(q,t)", kh, want);
        check_jones("left trefoil", kh, jones_chi_left_trefoil());

        /* Expect a Z_2 torsion summand (integral Khovanov of the trefoil) */
        int found_z2 = 0;
        for (int i = 0; i < tors.count; i++) {
            if (tors.items[i].d == 2)
                found_z2 = 1;
        }
        if (found_z2) {
            printf("  PASS  left trefoil has Z_2 torsion\n");
        } else {
            printf("  FAIL  left trefoil missing Z_2 torsion (count=%d)\n", tors.count);
            g_failures++;
        }

        bp_free(want);
        bp_free(kh);
        kh_torsion_report_free(&tors);
        knot_free(k);
        printf("\n");
    }

    /* ---- Right trefoil ---- */
    {
        printf("[right trefoil]\n");
        Knot *k = knot_create_right_trefoil();
        BivariatePoly *kh = kh_poincare(k);
        BivariatePoly *want = gold_right_trefoil();
        expect_poly("right trefoil Kh(q,t)", kh, want);
        check_jones("right trefoil", kh, jones_chi_right_trefoil());
        bp_free(want);
        bp_free(kh);
        knot_free(k);
        printf("\n");
    }

    /* ---- Mirror consistency: left should be (q,t)->(q^{-1},t^{-1}) of right ---- */
    {
        printf("[mirror symmetry]\n");
        Knot *L = knot_create_left_trefoil();
        Knot *R = knot_create_right_trefoil();
        BivariatePoly *khL = kh_poincare(L);
        BivariatePoly *khR = kh_poincare(R);
        BivariatePoly *mirrorR = bp_create();
        for (int i = 0; i < khR->num_terms; i++) {
            bp_add_term(mirrorR,
                        -khR->terms[i].q_exp,
                        -khR->terms[i].t_exp,
                        khR->terms[i].coeff);
        }
        expect_poly("left == mirror(right)", khL, mirrorR);
        bp_free(mirrorR);
        bp_free(khL);
        bp_free(khR);
        knot_free(L);
        knot_free(R);
        printf("\n");
    }

    if (g_failures == 0) {
        printf("=== ALL GOLD TESTS PASSED ===\n");
        return 0;
    }
    printf("=== %d FAILURE(S) ===\n", g_failures);
    return 1;
}
