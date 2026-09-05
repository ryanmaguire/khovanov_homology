#include "TopologyBridge.h"

#include <stdlib.h>

BivariatePoly *kh_poincare_from_rank_table(const KhRankTable *table)
{
    return kh_poincare_from_rank_table_with_torsion(table, NULL, NULL);
}

BivariatePoly *kh_poincare_from_rank_table_with_torsion(
    const KhRankTable *table,
    const KhTorsionTable *tors_table,
    KhTorsionReport *out_tors)
{
    BivariatePoly *poly = bp_create();
    if (!poly)
        return NULL;
    if (!table || table->h_max < table->h_min || table->q_max < table->q_min)
        return poly;
    if (!table->chain_dims)
        return poly;

    int h_span = table->h_max - table->h_min + 1;
    int q_span = table->q_max - table->q_min + 1;

    for (int hi = 0; hi < h_span; hi++) {
        if (!table->chain_dims[hi])
            continue;
        int h = table->h_min + hi;
        int true_r = h - table->n_minus;

        for (int qi = 0; qi < q_span; qi++) {
            int q = table->q_min + qi;
            int chain = table->chain_dims[hi][qi];
            int rank_out = 0;
            int rank_in = 0;

            if (table->image_ranks) {
                if (h < table->h_max && table->image_ranks[hi])
                    rank_out = table->image_ranks[hi][qi];
                if (hi > 0 && table->image_ranks[hi - 1])
                    rank_in = table->image_ranks[hi - 1][qi];
            }

            int beta = chain - rank_out - rank_in;
            if (beta > 0)
                bp_add_term(poly, q, true_r, beta);

            /* optional torsion on outgoing d_h, attributed to target degree */
            if (out_tors && tors_table && tors_table->tors_counts &&
                tors_table->tors_values && h < table->h_max &&
                tors_table->tors_counts[hi]) {
                int n = tors_table->tors_counts[hi][qi];
                int *vals = tors_table->tors_values[hi]
                                ? tors_table->tors_values[hi][qi]
                                : NULL;
                int target_r = (h + 1) - table->n_minus;
                for (int k = 0; k < n && vals; k++) {
                    if (vals[k] > 1) {
                        /* reuse KhPoincare torsion push via report API:
                         * expand capacity manually here to avoid depending
                         * on static helpers */
                        if (out_tors->count == out_tors->capacity) {
                            int nc = out_tors->capacity == 0 ? 8
                                                             : out_tors->capacity * 2;
                            KhTorsionSummand *ni = (KhTorsionSummand *)realloc(
                                out_tors->items,
                                (size_t)nc * sizeof(KhTorsionSummand));
                            if (!ni)
                                break;
                            out_tors->items = ni;
                            out_tors->capacity = nc;
                        }
                        out_tors->items[out_tors->count].r = target_r;
                        out_tors->items[out_tors->count].j = q;
                        out_tors->items[out_tors->count].d = vals[k];
                        out_tors->count++;
                    }
                }
            }
        }
    }

    return poly;
}
