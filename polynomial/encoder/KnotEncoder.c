#include "KnotEncoder.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/**
 * @brief Cross-platform popcount (number of 1-bits).
 *        Used for homological grading (h = pop_count(r)).
 */
int pop_count(unsigned int n)
{
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcount(n);
#elif defined(_MSC_VER)
    return __popcnt(n);
#else
    /* Portable fallback (Brian Kernighan algorithm) */
    int count = 0;
    while (n != 0) {
        n &= (n - 1);
        count++;
    }
    return count;
#endif
}

/* ========================================================================= */
/* Union-Find core logic for counting connected loops after smoothing        */
/* ========================================================================= */

/**
 * @brief Union-Find Find operation with path compression (highly efficient).
 */
static int uf_find(int* parent, int i) {
    int root = i;
    /* Find root */
    while (root != parent[root]) {
        root = parent[root];
    }
    /* Path compression: connect every node directly to root */
    int curr = i;
    while (curr != root) {
        int nxt = parent[curr];
        parent[curr] = root;
        curr = nxt;
    }
    return root;
}

/**
 * @brief Given a resolution state r, compute the number of independent loops
 *        formed after smoothing.
 *
 * @param knot  Knot containing Planar Diagram (PD) information
 * @param r     Current resolution bitmask
 * @return      Number of connected loops (circles)
 *
 * Note: Uses stack-allocated arrays to avoid malloc overhead for 2^m calls.
 *       Max edges = 4 * m (for m=19, only 76 edges needed).
 */
int knot_count_loops(Knot* knot, int r) {
    if (!knot || knot->m <= 0) return 0;

    /* 19-crossing knot has at most 76 edges; use fixed stack array */
    int max_edges = 4 * knot->m + 10;
    int parent[200];
    int active[200];

    /* Initialize Union-Find */
    for (int i = 0; i < max_edges; i++) {
        parent[i] = i;
        active[i] = 0;
    }

    /* Traverse all crossings and apply 0/1-smoothing */
    for (int k = 0; k < knot->m; k++) {
        int a = knot->pd[k][0], b = knot->pd[k][1];
        int c = knot->pd[k][2], d = knot->pd[k][3];

        /* Mark edges as active */
        active[a] = 1; active[b] = 1; active[c] = 1; active[d] = 1;

        int u1, v1, u2, v2;
        if ((r & (1 << k)) == 0) {
            /* 0-smoothing: connect (a,b) and (c,d) */
            u1 = a; v1 = b;
            u2 = c; v2 = d;
        } else {
            /* 1-smoothing: connect (a,d) and (b,c) */
            u1 = a; v1 = d;
            u2 = b; v2 = c;
        }

        /* Union(u1, v1) */
        int root_u1 = uf_find(parent, u1);
        int root_v1 = uf_find(parent, v1);
        if (root_u1 != root_v1) parent[root_v1] = root_u1;

        /* Union(u2, v2) */
        int root_u2 = uf_find(parent, u2);
        int root_v2 = uf_find(parent, v2);
        if (root_u2 != root_v2) parent[root_v2] = root_u2;
    }

    /* Count independent loops (active nodes that are their own root) */
    int num_loops = 0;
    for (int i = 0; i < max_edges; i++) {
        if (active[i] && parent[i] == i) {
            num_loops++;
        }
    }

    return num_loops;
}

/**
 * @brief Helper: compute writhe from PD (standard sign rule).
 */
static int compute_writhe(int pd[][4], int m) {
    int w = 0;
    for (int i = 0; i < m; i++) {
        if (pd[i][1] - pd[i][3] == 1 || pd[i][3] - pd[i][1] > 1) w++;
        else w--;
    }
    return w;
}

Knot* knot_from_pd(int pd[][4], int m) {
    Knot* k = malloc(sizeof(Knot));
    k->m = m;
    k->pd = malloc(m * sizeof(int*));
    for (int i = 0; i < m; i++) {
        k->pd[i] = malloc(4 * sizeof(int));
        for (int j = 0; j < 4; j++) k->pd[i][j] = pd[i][j];
    }
    k->writhe = compute_writhe(pd, m);
    return k;
}

void knot_free(Knot* k) {
    if (k) {
        for (int i = 0; i < k->m; i++) free(k->pd[i]);
        free(k->pd);
        free(k);
    }
}

/**
 * @brief Converts a resolution bitmask r into Arjun's SmoothingColumn structure.
 *        This is the key bridge between our encoder and Arjun's CobMatrix system.
 */
SmoothingColumn* resolution_to_smoothing_column(Knot* knot, int r) {
    int m = knot->m;
    SmoothingColumn* col = malloc(sizeof(SmoothingColumn));
    col->n = 2 * m; /* 2m edges */
    col->numbers = calloc(col->n, sizeof(int));
    col->smoothings = malloc(col->n * sizeof(Cap*));
    for (int k = 0; k < m; k++) {
        int a = knot->pd[k][0], b = knot->pd[k][1];
        int c = knot->pd[k][2], d = knot->pd[k][3];
        Cap* cap = Cap_create(2 * m, 0);
        if ((r & (1 << k)) == 0) { /* 0-smoothing */
            cap->pairings[a] = b; cap->pairings[b] = a;
            cap->pairings[c] = d; cap->pairings[d] = c;
        } else { /* 1-smoothing */
            cap->pairings[a] = d; cap->pairings[d] = a;
            cap->pairings[b] = c; cap->pairings[c] = b;
        }
        col->smoothings[k] = cap;
    }
    return col;
}

/**
 * @brief Constructs the CobMatrix representing the boundary operator ∂ between
 *        two SmoothingColumns. This is the core differential in the chain complex.
 */
CobMatrix* build_boundary_matrix(SmoothingColumn* source, SmoothingColumn* target, int r) {
    CobMatrix* m = CobMatrix_create(source, target, false);
    int m_cross = source->n / 2; /* number of crossings */
    for (int k = 0; k < m_cross; k++) {
        /* Create cobordism for this crossing */
        CannedCobordismImplData* impl = CannedCobordismImpl_create(
            source->smoothings[k], target->smoothings[k]);
        CannedCobordism* cc = CannedCobordismImpl_as_CannedCobordism(impl);
        /* Jordan-Wigner sign accumulation (paper 4.4) */
        int sign = 1;
        for (int j = 0; j < k; j++) {
            if (r & (1 << j)) sign = -sign;
        }
        LCCC* lc = LCCC_create(sign);
        /* Place in matrix (simplified diagonal for now; can be refined later) */
        CobMatrix_putEntry(m, k, k, lc);
    }
    return m;
}

BivariatePoly* compute_khovanov_polynomial(Knot* knot)
{
    /* Handle invalid or empty knot */
    if (!knot || knot->m == 0) {
        BivariatePoly* zero = bp_create();
        return zero;
    }

    int m = knot->m;
    BivariatePoly* poly = bp_create();

    /* Create full Komplex (the complete chain complex) */
    Komplex* komplex = Komplex_create(1 << m);

    /* Iterate over all adjacent resolutions and build the differentials */
    for (int r = 0; r < (1 << m) - 1; r++) {
        /* Generate SmoothingColumn for current resolution r */
        SmoothingColumn* source = resolution_to_smoothing_column(knot, r);
        /* Generate SmoothingColumn for next resolution (r+1) */
        SmoothingColumn* target = resolution_to_smoothing_column(knot, r + 1);

        /* Build the boundary matrix ∂ between these two smoothings */
        CobMatrix* boundary = build_boundary_matrix(source, target, r);

        /* Store the boundary matrix into Komplex */
        komplex->differentials[r] = boundary;

        /* Immediately attempt block reduction */
        Komplex_blockReductionLemma(komplex, r, 0, 0);

        /* Clean up temporary SmoothingColumn objects */
        for (int i = 0; i < source->n; i++) Cap_free(source->smoothings[i]);
        free(source->smoothings);
        free(source->numbers);
        free(source);

        for (int i = 0; i < target->n; i++) Cap_free(target->smoothings[i]);
        free(target->smoothings);
        free(target->numbers);
        free(target);
    }

    /* For now we compute the graded Euler characteristic as a fast approximation */
    for (int r = 0; r < (1 << m); r++) {
        int hom_deg = pop_count(r); /* homological degree i = |r| */

        /* Use the real loop count from Union-Find */
        int num_loops = knot_count_loops(knot, r);

        int quantum_deg = knot->writhe + hom_deg - 2 * num_loops;

        /* Simple sign for Euler characteristic contribution */
        int sign = (pop_count(r) % 2 == 0) ? 1 : -1;

        bp_add_term(poly, quantum_deg, hom_deg, sign);
    }

    /* Clean up Komplex */
    Komplex_free(komplex);
    return poly;
}
