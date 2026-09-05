#include "KnotEncoder.h"
#include "../complex/KhovanovKomplex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//debug helper
static void debug_report_knotencoder(const char *location, const char *msg, int a, int b, int c) {
    char command[2048];
    snprintf(command, sizeof(command),
             "sh -lc 'source .dbg/trefoil-poincare-hang.env 2>/dev/null || { DEBUG_SERVER_URL=\"http://127.0.0.1:7777/event\"; DEBUG_SESSION_ID=\"trefoil-poincare-hang\"; }; curl -sX POST \"$DEBUG_SERVER_URL\" -H \"Content-Type: application/json\" -d \"{\\\"sessionId\\\":\\\"$DEBUG_SESSION_ID\\\",\\\"runId\\\":\\\"pre-fix\\\",\\\"hypothesisId\\\":\\\"A\\\",\\\"location\\\":\\\"%s\\\",\\\"msg\\\":\\\"%s\\\",\\\"data\\\":{\\\"a\\\":%d,\\\"b\\\":%d,\\\"c\\\":%d}}\" >/dev/null 2>&1'",
             location, msg, a, b, c);
    (void)system(command);
}
//debug end

//popcount for homological grading

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
/* Union-Find core logic for counting connected loops after smoothing        */
/* ========================================================================= */

static int uf_find(int* parent, int i) {
    int root = i;
    while (root != parent[root]) {
        root = parent[root];
    }
    int curr = i;
    while (curr != root) {
        int nxt = parent[curr];
        parent[curr] = root;
        curr = nxt;
    }
    return root;
}

int knot_count_loops(Knot* knot, int r) {
    if (!knot || knot->m <= 0) return 0;

    int max_edges = 4 * knot->m + 10;
    int parent[max_edges];
    int active[max_edges];

    for (int i = 0; i < max_edges; i++) {
        parent[i] = i;
        active[i] = 0;
    }

    for (int k = 0; k < knot->m; k++) {
        int a = knot->pd[k][0], b = knot->pd[k][1];
        int c = knot->pd[k][2], d = knot->pd[k][3];

        active[a] = 1; active[b] = 1; active[c] = 1; active[d] = 1;

        int u1, v1, u2, v2;
        if ((r & (1 << k)) == 0) {
            u1 = a; v1 = b; u2 = c; v2 = d;
        } else {
            u1 = a; v1 = d; u2 = b; v2 = c;
        }

        int root_u1 = uf_find(parent, u1);
        int root_v1 = uf_find(parent, v1);
        if (root_u1 != root_v1) parent[root_v1] = root_u1;

        int root_u2 = uf_find(parent, u2);
        int root_v2 = uf_find(parent, v2);
        if (root_u2 != root_v2) parent[root_v2] = root_u2;
    }

    int num_loops = 0;
    for (int i = 0; i < max_edges; i++) {
        if (active[i] && parent[i] == i) {
            num_loops++;
        }
    }

    return num_loops;
}

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
        for (int j = 0; j < 4; j++) k->pd[i][j] = pd[i][j] - 1;
    }
    k->writhe = compute_writhe(pd, m);
    return k;
}

Knot* knot_create_trefoil(void) {
    int pd[3][4] = {
        {1, 4, 2, 5},
        {3, 6, 4, 1},
        {5, 2, 6, 3}
    };
    return knot_from_pd(pd, 3);
}

void knot_free(Knot* k) {
    if (k) {
        for (int i = 0; i < k->m; i++) free(k->pd[i]);
        free(k->pd);
        free(k);
    }
}

//computation pipeline
BivariatePoly* compute_khovanov_polynomial(Knot* knot) {
    if (!knot || knot->m == 0) return bp_create();

    BivariatePoly* poly = bp_create();
    
    //debug start
    debug_report_knotencoder("KnotEncoder.c:compute:start", "start", knot->m, knot->writhe, 0);
    //debug end
    KhovanovKomplex* kh_complex = KhovanovKomplex_alloc(knot);
    //debug alloc
    debug_report_knotencoder("KnotEncoder.c:alloc:done", kh_complex ? "alloc ok" : "alloc null", knot->m, knot->writhe, kh_complex != NULL);
    //debug end
    if (!kh_complex || !KhovanovKomplex_build(kh_complex)) {
        //build-failed
        debug_report_knotencoder("KnotEncoder.c:build:failed", "build fail", knot->m, knot->writhe, kh_complex != NULL);
        //debug end
        printf("Error: Failed to build KhovanovKomplex.\n");
        if (kh_complex) KhovanovKomplex_free(kh_complex);
        return poly;
    }
    //debug done
    debug_report_knotencoder("KnotEncoder.c:build:done", "build ok", knot->m, knot->writhe, 1);
    //debug end

    //get graded euler
    for (int h = 0; h <= knot->m; h++) {
        SmoothingColumn* col = kh_complex->columns[h];
        int sign = (h % 2 == 0) ? 1 : -1;

        for (int i = 0; i < col->n; i++) {
            int num_loops = col->numbers[i];
            int quantum_deg = knot->writhe + h - 2 * num_loops;
            bp_add_term(poly, quantum_deg, h, sign);
        }
    }

    //clean
    KhovanovKomplex_free(kh_complex);
    return poly;
}
