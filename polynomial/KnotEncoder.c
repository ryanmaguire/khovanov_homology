#include "KnotEncoder.h"
#include "../homology/KhPoincare.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* popcount for homological grading of a resolution bitmask */
int pop_count(unsigned int n)
{
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcount(n);
#elif defined(_MSC_VER)
    return __popcnt(n);
#else
    int count = 0;
    while (n != 0) {
        n &= (n - 1);
        count++;
    }
    return count;
#endif
}

/* ========================================================================= */
/* Union-Find: count connected loops after smoothing                         */
/* ========================================================================= */

static int uf_find(int *parent, int i)
{
    int root = i;
    while (root != parent[root])
        root = parent[root];
    int curr = i;
    while (curr != root) {
        int nxt = parent[curr];
        parent[curr] = root;
        curr = nxt;
    }
    return root;
}

int knot_count_loops(Knot *knot, int r)
{
    if (!knot || knot->m <= 0)
        return 0;

    int max_edges = 4 * knot->m + 10;
    int *parent = (int *)malloc((size_t)max_edges * sizeof(int));
    int *active = (int *)malloc((size_t)max_edges * sizeof(int));
    if (!parent || !active) {
        free(parent);
        free(active);
        return 0;
    }

    for (int i = 0; i < max_edges; i++) {
        parent[i] = i;
        active[i] = 0;
    }

    for (int k = 0; k < knot->m; k++) {
        int a = knot->pd[k][0], b = knot->pd[k][1];
        int c = knot->pd[k][2], d = knot->pd[k][3];

        if (a >= 0 && a < max_edges) active[a] = 1;
        if (b >= 0 && b < max_edges) active[b] = 1;
        if (c >= 0 && c < max_edges) active[c] = 1;
        if (d >= 0 && d < max_edges) active[d] = 1;

        int u1, v1, u2, v2;
        if ((r & (1 << k)) == 0) {
            u1 = a;
            v1 = b;
            u2 = c;
            v2 = d;
        } else {
            u1 = a;
            v1 = d;
            u2 = b;
            v2 = c;
        }

        int root_u1 = uf_find(parent, u1);
        int root_v1 = uf_find(parent, v1);
        if (root_u1 != root_v1)
            parent[root_v1] = root_u1;

        int root_u2 = uf_find(parent, u2);
        int root_v2 = uf_find(parent, v2);
        if (root_u2 != root_v2)
            parent[root_v2] = root_u2;
    }

    int num_loops = 0;
    for (int i = 0; i < max_edges; i++) {
        if (active[i] && parent[i] == i)
            num_loops++;
    }

    free(parent);
    free(active);
    return num_loops;
}

static int compute_writhe(int pd[][4], int m)
{
    int w = 0;
    for (int i = 0; i < m; i++) {
        if (pd[i][1] - pd[i][3] == 1 || pd[i][3] - pd[i][1] > 1)
            w++;
        else
            w--;
    }
    return w;
}

Knot *knot_from_pd(int pd[][4], int m)
{
    Knot *k = (Knot *)malloc(sizeof(Knot));
    if (!k)
        return NULL;
    k->m = m;
    k->edges = NULL;
    k->pd = (int **)malloc((size_t)m * sizeof(int *));
    if (!k->pd) {
        free(k);
        return NULL;
    }
    for (int i = 0; i < m; i++) {
        k->pd[i] = (int *)malloc(4 * sizeof(int));
        if (!k->pd[i]) {
            for (int j = 0; j < i; j++)
                free(k->pd[j]);
            free(k->pd);
            free(k);
            return NULL;
        }
        /* convert 1-based PD labels to 0-based */
        for (int j = 0; j < 4; j++)
            k->pd[i][j] = pd[i][j] - 1;
    }
    k->writhe = compute_writhe(pd, m);
    return k;
}

Knot *knot_create_trefoil(void)
{
    int pd[3][4] = {
        {1, 4, 2, 5},
        {3, 6, 4, 1},
        {5, 2, 6, 3}
    };
    return knot_from_pd(pd, 3);
}

void knot_free(Knot *k)
{
    if (!k)
        return;
    if (k->pd) {
        for (int i = 0; i < k->m; i++)
            free(k->pd[i]);
        free(k->pd);
    }
    free(k->edges);
    free(k);
}

/* ========================================================================= */
/* Production Poincaré polynomial — true homology                            */
/* ========================================================================= */

BivariatePoly *compute_khovanov_polynomial(Knot *knot)
{
    /* Delegate to enhanced-state + Smith engine (step 3). */
    return kh_poincare(knot);
}

/* ========================================================================= */
/* Legacy: signed state sum over the cube (graded Euler, NOT homology)       */
/* ========================================================================= */

BivariatePoly *compute_graded_euler_state_sum(Knot *knot)
{
    if (!knot || knot->m == 0)
        return bp_create();

    BivariatePoly *poly = bp_create();
    if (!poly)
        return NULL;

    int m = knot->m;
    int num_states = 1 << m;

    for (int r = 0; r < num_states; r++) {
        int h = pop_count((unsigned int)r);
        int num_loops = knot_count_loops(knot, r);
        int sign = (h % 2 == 0) ? 1 : -1;
        /* cube quantum degree (unnormalized): w + h - 2ℓ */
        int quantum_deg = knot->writhe + h - 2 * num_loops;
        bp_add_term(poly, quantum_deg, h, sign);
    }

    return poly;
}
