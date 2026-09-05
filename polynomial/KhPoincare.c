#include "KhPoincare.h"
#include "../../IntegerMatrix.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- internal types (from trefoil_integer_poincare) ---- */

typedef struct {
    int count;
    int capacity;
    int *resolutions;
    unsigned int *masks;
} BasisBucket;

typedef struct {
    int n_occ;
    int n_circles;
    int *comp_of_occ;
} ResolutionInfo;

/* ---- small helpers ---- */

static int pop_count_u32(uint32_t x)
{
    int c = 0;
    while (x) {
        x &= (x - 1);
        c++;
    }
    return c;
}

static int jordan_wigner_sign(int src_res, int dst_res)
{
    int diff = src_res ^ dst_res;
    if (diff == 0 || (diff & (diff - 1)) != 0 || (diff & dst_res) == 0)
        return 0;
    int mask = diff - 1;
    int prior_ones = src_res & mask;
    return (pop_count_u32((uint32_t)prior_ones) % 2 == 0) ? 1 : -1;
}

static int dsu_find(int *parent, int x)
{
    if (parent[x] != x)
        parent[x] = dsu_find(parent, parent[x]);
    return parent[x];
}

static void dsu_union(int *parent, int a, int b)
{
    int ra = dsu_find(parent, a);
    int rb = dsu_find(parent, b);
    if (ra != rb)
        parent[rb] = ra;
}

static void bucket_push(BasisBucket *bucket, int resolution, unsigned int mask)
{
    if (bucket->count == bucket->capacity) {
        int new_capacity = bucket->capacity == 0 ? 8 : bucket->capacity * 2;
        int *new_resolutions =
            (int *)realloc(bucket->resolutions, (size_t)new_capacity * sizeof(int));
        unsigned int *new_masks =
            (unsigned int *)realloc(bucket->masks, (size_t)new_capacity * sizeof(unsigned int));
        if (!new_resolutions || !new_masks) {
            fprintf(stderr, "KhPoincare: out of memory growing basis bucket\n");
            free(new_resolutions);
            free(new_masks);
            return;
        }
        bucket->resolutions = new_resolutions;
        bucket->masks = new_masks;
        bucket->capacity = new_capacity;
    }
    bucket->resolutions[bucket->count] = resolution;
    bucket->masks[bucket->count] = mask;
    bucket->count++;
}

static int bucket_find(const BasisBucket *bucket, int resolution, unsigned int mask)
{
    for (int i = 0; i < bucket->count; i++) {
        if (bucket->resolutions[i] == resolution && bucket->masks[i] == mask)
            return i;
    }
    return -1;
}

static void bucket_free(BasisBucket *bucket)
{
    free(bucket->resolutions);
    free(bucket->masks);
    bucket->resolutions = NULL;
    bucket->masks = NULL;
    bucket->count = 0;
    bucket->capacity = 0;
}

static void compute_label_occurrences(const Knot *knot, int *first_occ, int *second_occ)
{
    int labels = 2 * knot->m;
    for (int i = 0; i < labels; i++) {
        first_occ[i] = -1;
        second_occ[i] = -1;
    }
    for (int c = 0; c < knot->m; c++) {
        for (int p = 0; p < 4; p++) {
            int label = knot->pd[c][p];
            int occ = 4 * c + p;
            if (label < 0 || label >= labels)
                continue;
            if (first_occ[label] < 0)
                first_occ[label] = occ;
            else
                second_occ[label] = occ;
        }
    }
}

static void compute_resolution_info(const Knot *knot,
                                    const int *first_occ,
                                    const int *second_occ,
                                    int resolution,
                                    ResolutionInfo *info)
{
    int n_occ = 4 * knot->m;
    int labels = 2 * knot->m;
    int *parent = (int *)malloc((size_t)n_occ * sizeof(int));
    int *root_to_comp = (int *)malloc((size_t)n_occ * sizeof(int));
    if (!parent || !root_to_comp) {
        free(parent);
        free(root_to_comp);
        info->n_occ = 0;
        info->n_circles = 0;
        info->comp_of_occ = NULL;
        return;
    }

    for (int i = 0; i < n_occ; i++) {
        parent[i] = i;
        root_to_comp[i] = -1;
    }

    for (int label = 0; label < labels; label++) {
        if (first_occ[label] >= 0 && second_occ[label] >= 0)
            dsu_union(parent, first_occ[label], second_occ[label]);
    }

    for (int c = 0; c < knot->m; c++) {
        int o0 = 4 * c + 0, o1 = 4 * c + 1, o2 = 4 * c + 2, o3 = 4 * c + 3;
        if ((resolution & (1 << c)) == 0) {
            dsu_union(parent, o0, o1);
            dsu_union(parent, o2, o3);
        } else {
            dsu_union(parent, o0, o3);
            dsu_union(parent, o1, o2);
        }
    }

    info->n_occ = n_occ;
    info->n_circles = 0;
    info->comp_of_occ = (int *)malloc((size_t)n_occ * sizeof(int));
    if (!info->comp_of_occ) {
        free(parent);
        free(root_to_comp);
        info->n_circles = 0;
        return;
    }

    for (int occ = 0; occ < n_occ; occ++) {
        int root = dsu_find(parent, occ);
        if (root_to_comp[root] < 0)
            root_to_comp[root] = info->n_circles++;
        info->comp_of_occ[occ] = root_to_comp[root];
    }

    free(parent);
    free(root_to_comp);
}

static void free_resolution_info(ResolutionInfo *info)
{
    free(info->comp_of_occ);
    info->comp_of_occ = NULL;
}

static void add_matrix_entry(Mat *differential,
                             const BasisBucket *target_bucket,
                             int target_resolution,
                             unsigned int target_mask,
                             int target_coeff,
                             int source_col)
{
    int row = bucket_find(target_bucket, target_resolution, target_mask);
    if (row >= 0)
        differential->matrix[row][source_col] += target_coeff;
}

static void emit_edge_terms(const ResolutionInfo *src,
                            const ResolutionInfo *dst,
                            int target_resolution,
                            unsigned int source_mask,
                            int cube_sign,
                            const BasisBucket *target_bucket,
                            Mat *differential,
                            int source_col)
{
    int src_k = src->n_circles;
    int dst_k = dst->n_circles;
    if (src_k <= 0 || dst_k <= 0)
        return;

    int *overlap_counts_src =
        (int *)calloc((size_t)src_k * (size_t)dst_k, sizeof(int));
    int *dst_from_src = (int *)malloc((size_t)dst_k * sizeof(int));
    int *src_from_dst = (int *)malloc((size_t)src_k * sizeof(int));
    if (!overlap_counts_src || !dst_from_src || !src_from_dst) {
        free(overlap_counts_src);
        free(dst_from_src);
        free(src_from_dst);
        return;
    }

    for (int i = 0; i < dst_k; i++)
        dst_from_src[i] = -1;
    for (int i = 0; i < src_k; i++)
        src_from_dst[i] = -1;
    for (int occ = 0; occ < src->n_occ; occ++)
        overlap_counts_src[src->comp_of_occ[occ] * dst_k + dst->comp_of_occ[occ]] = 1;

    if (dst_k == src_k - 1) {
        /* merge: two source circles → one target circle */
        int merge_dst = -1, merge_src_a = -1, merge_src_b = -1;
        for (int d = 0; d < dst_k; d++) {
            int count = 0, first = -1, second = -1;
            for (int s = 0; s < src_k; s++) {
                if (overlap_counts_src[s * dst_k + d]) {
                    if (count == 0)
                        first = s;
                    else if (count == 1)
                        second = s;
                    count++;
                }
            }
            if (count == 2) {
                merge_dst = d;
                merge_src_a = first;
                merge_src_b = second;
            } else if (count == 1) {
                dst_from_src[d] = first;
            }
        }

        if (merge_dst < 0 || merge_src_a < 0 || merge_src_b < 0)
            goto cleanup;

        int bit_a = (int)((source_mask >> merge_src_a) & 1u);
        int bit_b = (int)((source_mask >> merge_src_b) & 1u);
        if (bit_a == 1 && bit_b == 1)
            goto cleanup; /* m(X,X) = 0 */

        int merged_bit = (bit_a == 0 && bit_b == 0) ? 0 : 1;
        unsigned int target_mask = 0u;
        for (int d = 0; d < dst_k; d++) {
            int bit = (d == merge_dst)
                          ? merged_bit
                          : (int)((source_mask >> dst_from_src[d]) & 1u);
            if (bit)
                target_mask |= (1u << d);
        }
        add_matrix_entry(differential, target_bucket, target_resolution,
                         target_mask, cube_sign, source_col);

    } else if (dst_k == src_k + 1) {
        /* split: one source circle → two target circles */
        int split_src = -1, split_dst_a = -1, split_dst_b = -1;
        for (int s = 0; s < src_k; s++) {
            int count = 0, first = -1, second = -1;
            for (int d = 0; d < dst_k; d++) {
                if (overlap_counts_src[s * dst_k + d]) {
                    if (count == 0)
                        first = d;
                    else if (count == 1)
                        second = d;
                    count++;
                }
            }
            if (count == 2) {
                split_src = s;
                split_dst_a = first;
                split_dst_b = second;
            } else if (count == 1) {
                src_from_dst[s] = first;
            }
        }

        if (split_src < 0 || split_dst_a < 0 || split_dst_b < 0)
            goto cleanup;

        int split_bit = (int)((source_mask >> split_src) & 1u);
        int term_count = (split_bit == 0) ? 2 : 1;
        int split_terms[2][2];
        if (split_bit == 0) {
            /* Δ(1) = 1⊗X + X⊗1 */
            split_terms[0][0] = 0;
            split_terms[0][1] = 1;
            split_terms[1][0] = 1;
            split_terms[1][1] = 0;
        } else {
            /* Δ(X) = X⊗X */
            split_terms[0][0] = 1;
            split_terms[0][1] = 1;
        }

        for (int term = 0; term < term_count; term++) {
            unsigned int target_mask = 0u;
            for (int d = 0; d < dst_k; d++) {
                int bit = 0;
                if (d == split_dst_a)
                    bit = split_terms[term][0];
                else if (d == split_dst_b)
                    bit = split_terms[term][1];
                else {
                    int unique_src = -1;
                    for (int s = 0; s < src_k; s++) {
                        if (overlap_counts_src[s * dst_k + d]) {
                            unique_src = s;
                            break;
                        }
                    }
                    if (unique_src >= 0)
                        bit = (int)((source_mask >> unique_src) & 1u);
                }
                if (bit)
                    target_mask |= (1u << d);
            }
            add_matrix_entry(differential, target_bucket, target_resolution,
                             target_mask, cube_sign, source_col);
        }
    }

cleanup:
    free(overlap_counts_src);
    free(dst_from_src);
    free(src_from_dst);
}

/* ---- torsion report ---- */

void kh_torsion_report_init(KhTorsionReport *rep)
{
    if (!rep)
        return;
    rep->items = NULL;
    rep->count = 0;
    rep->capacity = 0;
}

void kh_torsion_report_free(KhTorsionReport *rep)
{
    if (!rep)
        return;
    free(rep->items);
    rep->items = NULL;
    rep->count = 0;
    rep->capacity = 0;
}

static void torsion_push(KhTorsionReport *rep, int r, int j, int d)
{
    if (!rep || d <= 1)
        return;
    if (rep->count == rep->capacity) {
        int nc = rep->capacity == 0 ? 8 : rep->capacity * 2;
        KhTorsionSummand *ni =
            (KhTorsionSummand *)realloc(rep->items, (size_t)nc * sizeof(KhTorsionSummand));
        if (!ni)
            return;
        rep->items = ni;
        rep->capacity = nc;
    }
    rep->items[rep->count].r = r;
    rep->items[rep->count].j = j;
    rep->items[rep->count].d = d;
    rep->count++;
}

/* ---- main computation ---- */

BivariatePoly *kh_poincare_with_torsion(const Knot *knot, KhTorsionReport *tors)
{
    BivariatePoly *empty = bp_create();
    if (!knot || knot->m < 0)
        return empty;
    if (knot->m == 0) {
        /* unknot-like empty diagram: return empty; caller may special-case */
        return empty;
    }
    if (knot->m > 16) {
        fprintf(stderr, "KhPoincare: m=%d exceeds full-cube limit (16)\n", knot->m);
        return empty;
    }

    int m = knot->m;
    int num_states = 1 << m;
    int n_plus = (m + knot->writhe) / 2;
    int n_minus = (m - knot->writhe) / 2;

    int *first_occ = (int *)malloc((size_t)(2 * m) * sizeof(int));
    int *second_occ = (int *)malloc((size_t)(2 * m) * sizeof(int));
    ResolutionInfo *states =
        (ResolutionInfo *)calloc((size_t)num_states, sizeof(ResolutionInfo));
    if (!first_occ || !second_occ || !states) {
        free(first_occ);
        free(second_occ);
        free(states);
        return empty;
    }

    compute_label_occurrences(knot, first_occ, second_occ);
    for (int r = 0; r < num_states; r++)
        compute_resolution_info(knot, first_occ, second_occ, r, &states[r]);

    int min_q = INT_MAX, max_q = INT_MIN;
    for (int r = 0; r < num_states; r++) {
        int h = pop_count_u32((uint32_t)r);
        int circles = states[r].n_circles;
        if (circles < 0 || circles > 20)
            continue;
        unsigned int basis_count = 1u << circles;
        for (unsigned int mask = 0; mask < basis_count; mask++) {
            int q = h + circles - 2 * pop_count_u32(mask) + n_plus - 2 * n_minus;
            if (q < min_q)
                min_q = q;
            if (q > max_q)
                max_q = q;
        }
    }

    if (min_q > max_q) {
        for (int r = 0; r < num_states; r++)
            free_resolution_info(&states[r]);
        free(states);
        free(first_occ);
        free(second_occ);
        return empty;
    }

    int q_count = max_q - min_q + 1;
    BasisBucket **buckets =
        (BasisBucket **)malloc((size_t)(m + 1) * sizeof(BasisBucket *));
    int **ranks = (int **)malloc((size_t)m * sizeof(int *));
    int **betti = (int **)malloc((size_t)(m + 1) * sizeof(int *));

    if (!buckets || !ranks || !betti) {
        free(buckets);
        free(ranks);
        free(betti);
        for (int r = 0; r < num_states; r++)
            free_resolution_info(&states[r]);
        free(states);
        free(first_occ);
        free(second_occ);
        return empty;
    }

    for (int h = 0; h <= m; h++) {
        buckets[h] = (BasisBucket *)calloc((size_t)q_count, sizeof(BasisBucket));
        betti[h] = (int *)calloc((size_t)q_count, sizeof(int));
        if (!buckets[h] || !betti[h]) {
            /* partial cleanup below via same loops */
        }
    }
    for (int h = 0; h < m; h++)
        ranks[h] = (int *)calloc((size_t)q_count, sizeof(int));

    for (int r = 0; r < num_states; r++) {
        int h = pop_count_u32((uint32_t)r);
        int circles = states[r].n_circles;
        if (circles < 0 || circles > 20)
            continue;
        unsigned int basis_count = 1u << circles;
        for (unsigned int mask = 0; mask < basis_count; mask++) {
            int q = h + circles - 2 * pop_count_u32(mask) + n_plus - 2 * n_minus;
            bucket_push(&buckets[h][q - min_q], r, mask);
        }
    }

    /* differentials d_h : C_h → C_{h+1}, quantum-preserving */
    for (int h = 0; h < m; h++) {
        for (int q = min_q; q <= max_q; q++) {
            BasisBucket *source_bucket = &buckets[h][q - min_q];
            BasisBucket *target_bucket = &buckets[h + 1][q - min_q];
            if (source_bucket->count == 0) {
                ranks[h][q - min_q] = 0;
                continue;
            }

            Mat *differential =
                createMat(target_bucket->count, source_bucket->count);
            if (!differential) {
                ranks[h][q - min_q] = 0;
                continue;
            }

            for (int col = 0; col < source_bucket->count; col++) {
                int src_resolution = source_bucket->resolutions[col];
                unsigned int src_mask = source_bucket->masks[col];
                for (int crossing = 0; crossing < m; crossing++) {
                    if ((src_resolution & (1 << crossing)) != 0)
                        continue;
                    int dst_resolution = src_resolution | (1 << crossing);
                    int sign = jordan_wigner_sign(src_resolution, dst_resolution);
                    if (sign == 0)
                        continue;
                    emit_edge_terms(&states[src_resolution],
                                    &states[dst_resolution], dst_resolution,
                                    src_mask, sign, target_bucket, differential,
                                    col);
                }
            }

            int free_rank = 0;
            int tors_buf[32];
            int n_tors = 0;
            toSmithForm_extract(differential, &free_rank, tors_buf, 32, &n_tors);
            ranks[h][q - min_q] = free_rank;

            /* Torsion attributed to target degree (h+1, q), normalized r */
            if (tors) {
                int true_r = (h + 1) - n_minus;
                for (int i = 0; i < n_tors; i++)
                    torsion_push(tors, true_r, q, tors_buf[i]);
            }

            freeMat(differential);
        }
    }

    /* Betti: β_h = dim C_h − rank(d_h) − rank(d_{h−1}) */
    for (int h = 0; h <= m; h++) {
        for (int q = min_q; q <= max_q; q++) {
            int chain_rank = buckets[h][q - min_q].count;
            int rank_out = (h < m) ? ranks[h][q - min_q] : 0;
            int rank_in = (h > 0) ? ranks[h - 1][q - min_q] : 0;
            int b = chain_rank - rank_out - rank_in;
            betti[h][q - min_q] = (b > 0) ? b : 0;
        }
    }

    /* Assemble polynomial in normalized (r, j) = (h − n_-, q) */
    BivariatePoly *poly = bp_create();
    if (poly) {
        for (int h = 0; h <= m; h++) {
            int true_r = h - n_minus;
            for (int q = min_q; q <= max_q; q++) {
                int rank = betti[h][q - min_q];
                if (rank > 0)
                    bp_add_term(poly, q, true_r, rank);
            }
        }
    } else {
        poly = empty;
        empty = NULL;
    }

    /* cleanup */
    for (int r = 0; r < num_states; r++)
        free_resolution_info(&states[r]);
    free(states);
    free(first_occ);
    free(second_occ);

    for (int h = 0; h <= m; h++) {
        if (buckets[h]) {
            for (int qi = 0; qi < q_count; qi++)
                bucket_free(&buckets[h][qi]);
            free(buckets[h]);
        }
        free(betti[h]);
    }
    free(buckets);
    free(betti);
    for (int h = 0; h < m; h++)
        free(ranks[h]);
    free(ranks);

    if (empty)
        bp_free(empty);
    return poly ? poly : bp_create();
}

BivariatePoly *kh_poincare(const Knot *knot)
{
    return kh_poincare_with_torsion(knot, NULL);
}

Knot *knot_create_left_trefoil(void)
{
    /* 0-based PD labels, writhe -3 (matches trefoil_integer_poincare demo) */
    static const int pd_data[3][4] = {
        {4, 1, 3, 0},
        {2, 5, 1, 4},
        {0, 3, 5, 2}
    };
    Knot *knot = (Knot *)malloc(sizeof(Knot));
    if (!knot)
        return NULL;
    knot->m = 3;
    knot->writhe = -3;
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
