#include "../encoder/KnotEncoder.h"
#include "../../IntegerMatrix.h"
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

static int pop_count_u32(uint32_t x) {
    int c = 0;
    while (x) { x &= (x - 1); c++; }
    return c;
}

static int jordan_wigner_sign(int src_res, int dst_res) {
    int diff = src_res ^ dst_res;
    if (diff == 0 || (diff & (diff - 1)) != 0 || (diff & dst_res) == 0) return 0;
    int mask = diff - 1;
    int prior_ones = src_res & mask;
    return (pop_count_u32((uint32_t)prior_ones) % 2 == 0) ? 1 : -1;
}

static int dsu_find(int *parent, int x) {
    if (parent[x] != x) parent[x] = dsu_find(parent, parent[x]);
    return parent[x];
}

static void dsu_union(int *parent, int a, int b) {
    int ra = dsu_find(parent, a);
    int rb = dsu_find(parent, b);
    if (ra != rb) parent[rb] = ra;
}

static void bucket_push(BasisBucket *bucket, int resolution, unsigned int mask) {
    if (bucket->count == bucket->capacity) {
        int new_capacity = bucket->capacity == 0 ? 8 : bucket->capacity * 2;
        int *new_resolutions = (int *)realloc(bucket->resolutions, (size_t)new_capacity * sizeof(int));
        unsigned int *new_masks = (unsigned int *)realloc(bucket->masks, (size_t)new_capacity * sizeof(unsigned int));
        if (!new_resolutions || !new_masks) {
            fprintf(stderr, "Out of memory while growing bucket.\n");
            exit(2);
        }
        bucket->resolutions = new_resolutions;
        bucket->masks = new_masks;
        bucket->capacity = new_capacity;
    }
    bucket->resolutions[bucket->count] = resolution;
    bucket->masks[bucket->count] = mask;
    bucket->count++;
}

static int bucket_find(const BasisBucket *bucket, int resolution, unsigned int mask) {
    for (int i = 0; i < bucket->count; i++) {
        if (bucket->resolutions[i] == resolution && bucket->masks[i] == mask) return i;
    }
    return -1;
}

static void bucket_free(BasisBucket *bucket) {
    free(bucket->resolutions);
    free(bucket->masks);
}

static void compute_label_occurrences(Knot *knot, int *first_occ, int *second_occ) {
    int labels = 2 * knot->m;
    for (int i = 0; i < labels; i++) {
        first_occ[i] = -1;
        second_occ[i] = -1;
    }
    for (int c = 0; c < knot->m; c++) {
        for (int p = 0; p < 4; p++) {
            int label = knot->pd[c][p];
            int occ = 4 * c + p;
            if (first_occ[label] < 0) first_occ[label] = occ;
            else second_occ[label] = occ;
        }
    }
}

static void compute_resolution_info(Knot *knot, const int *first_occ, const int *second_occ, int resolution, ResolutionInfo *info) {
    int n_occ = 4 * knot->m;
    int labels = 2 * knot->m;
    int *parent = (int *)malloc((size_t)n_occ * sizeof(int));
    int *root_to_comp = (int *)malloc((size_t)n_occ * sizeof(int));
    
    for (int i = 0; i < n_occ; i++) {
        parent[i] = i;
        root_to_comp[i] = -1;
    }

    for (int label = 0; label < labels; label++) dsu_union(parent, first_occ[label], second_occ[label]);

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

    for (int occ = 0; occ < n_occ; occ++) {
        int root = dsu_find(parent, occ);
        if (root_to_comp[root] < 0) root_to_comp[root] = info->n_circles++;
        info->comp_of_occ[occ] = root_to_comp[root];
    }

    free(parent);
    free(root_to_comp);
}

static void free_resolution_info(ResolutionInfo *info) {
    free(info->comp_of_occ);
    info->comp_of_occ = NULL;
}

static void add_matrix_entry(Mat *differential, const BasisBucket *target_bucket, int target_resolution, unsigned int target_mask, int target_coeff, int source_col) {
    int row = bucket_find(target_bucket, target_resolution, target_mask);
    if (row >= 0) differential->matrix[row][source_col] += target_coeff;
}

static void emit_edge_terms(const ResolutionInfo *src, const ResolutionInfo *dst, int target_resolution, unsigned int source_mask, int cube_sign, const BasisBucket *target_bucket, Mat *differential, int source_col) {
    int src_k = src->n_circles;
    int dst_k = dst->n_circles;
    int *overlap_counts_src = (int *)calloc((size_t)src_k * (size_t)dst_k, sizeof(int));
    int *dst_from_src = (int *)malloc((size_t)dst_k * sizeof(int));
    int *src_from_dst = (int *)malloc((size_t)src_k * sizeof(int));

    for (int i = 0; i < dst_k; i++) dst_from_src[i] = -1;
    for (int i = 0; i < src_k; i++) src_from_dst[i] = -1;
    for (int occ = 0; occ < src->n_occ; occ++) overlap_counts_src[src->comp_of_occ[occ] * dst_k + dst->comp_of_occ[occ]] = 1;

    if (dst_k == src_k - 1) {
        int merge_dst = -1, merge_src_a = -1, merge_src_b = -1;
        for (int d = 0; d < dst_k; d++) {
            int count = 0, first = -1, second = -1;
            for (int s = 0; s < src_k; s++) {
                if (overlap_counts_src[s * dst_k + d]) {
                    if (count == 0) first = s; else if (count == 1) second = s;
                    count++;
                }
            }
            if (count == 2) { merge_dst = d; merge_src_a = first; merge_src_b = second; }
            else if (count == 1) dst_from_src[d] = first;
        }

        int bit_a = (int)((source_mask >> merge_src_a) & 1u);
        int bit_b = (int)((source_mask >> merge_src_b) & 1u);
        if (bit_a == 1 && bit_b == 1) goto cleanup;

        int merged_bit = (bit_a == 0 && bit_b == 0) ? 0 : 1;
        unsigned int target_mask = 0u;
        for (int d = 0; d < dst_k; d++) {
            int bit = (d == merge_dst) ? merged_bit : (int)((source_mask >> dst_from_src[d]) & 1u);
            if (bit) target_mask |= (1u << d);
        }
        add_matrix_entry(differential, target_bucket, target_resolution, target_mask, cube_sign, source_col);

    } else if (dst_k == src_k + 1) {
        int split_src = -1, split_dst_a = -1, split_dst_b = -1;
        for (int s = 0; s < src_k; s++) {
            int count = 0, first = -1, second = -1;
            for (int d = 0; d < dst_k; d++) {
                if (overlap_counts_src[s * dst_k + d]) {
                    if (count == 0) first = d; else if (count == 1) second = d;
                    count++;
                }
            }
            if (count == 2) { split_src = s; split_dst_a = first; split_dst_b = second; }
            else if (count == 1) src_from_dst[s] = first;
        }

        int split_bit = (int)((source_mask >> split_src) & 1u);
        int term_count = split_bit == 0 ? 2 : 1;
        int split_terms[2][2];
        if (split_bit == 0) { split_terms[0][0] = 0; split_terms[0][1] = 1; split_terms[1][0] = 1; split_terms[1][1] = 0; }
        else { split_terms[0][0] = 1; split_terms[0][1] = 1; }

        for (int term = 0; term < term_count; term++) {
            unsigned int target_mask = 0u;
            for (int d = 0; d < dst_k; d++) {
                int bit = 0;
                if (d == split_dst_a) bit = split_terms[term][0];
                else if (d == split_dst_b) bit = split_terms[term][1];
                else {
                    int unique_src = -1;
                    for (int s = 0; s < src_k; s++) if (overlap_counts_src[s * dst_k + d]) { unique_src = s; break; }
                    bit = (int)((source_mask >> unique_src) & 1u);
                }
                if (bit) target_mask |= (1u << d);
            }
            add_matrix_entry(differential, target_bucket, target_resolution, target_mask, cube_sign, source_col);
        }
    }

cleanup:
    free(overlap_counts_src);
    free(dst_from_src);
    free(src_from_dst);
}

static void print_bigraded_poincare(int **betti, int m, int min_q, int max_q, int n_minus, int *torsion_counts, int *torsion_values, int q_count) {
    bool first = true;
    printf("P(q,t) = ");
    for (int h = 0; h <= m; h++) {
        for (int q = min_q; q <= max_q; q++) {
            int idx = h * q_count + (q - min_q);
            int free_rank = betti[h][q - min_q];
            int t_count = torsion_counts[idx];
            
            if (free_rank == 0 && t_count == 0) continue;
            int true_h = h - n_minus;
            
            if (free_rank > 0) {
                if (!first) printf(" + ");
                first = false;
                if (free_rank != 1) printf("%d", free_rank);
                if (q == 0) printf("q^0"); else printf("q^%d", q);
                if (true_h != 0) printf("t^%d", true_h);
            }
            
            for (int i = 0; i < t_count; i++) {
                if (!first) printf(" + ");
                first = false;
                printf("Z_%d", torsion_values[idx * 10 + i]);
                if (q == 0) printf("q^0"); else printf("q^%d", q);
                if (true_h != 0) printf("t^%d", true_h);
            }
        }
    }
    if (first) printf("0");
    printf("\n");
}

static Knot *create_left_trefoil_knot(void) {
    static const int pd_data[3][4] = { {4, 1, 3, 0}, {2, 5, 1, 4}, {0, 3, 5, 2} };
    Knot *knot = (Knot *)malloc(sizeof(Knot));
    knot->m = 3; knot->writhe = -3; knot->edges = NULL;
    knot->pd = (int **)malloc((size_t)knot->m * sizeof(int *));
    for (int i = 0; i < knot->m; i++) {
        knot->pd[i] = (int *)malloc(4 * sizeof(int));
        for (int j = 0; j < 4; j++) knot->pd[i][j] = pd_data[i][j];
    }
    return knot;
}

static void free_local_knot(Knot *knot) {
    if (!knot) return;
    for (int i = 0; i < knot->m; i++) free(knot->pd[i]);
    free(knot->pd); free(knot->edges); free(knot);
}

int main(void) {
    Knot *knot = create_left_trefoil_knot();
    int m = knot->m;
    int num_states = 1 << m;
    int n_plus = (knot->m + knot->writhe) / 2;
    int n_minus = (knot->m - knot->writhe) / 2;

    int *first_occ = (int *)malloc((size_t)(2 * m) * sizeof(int));
    int *second_occ = (int *)malloc((size_t)(2 * m) * sizeof(int));
    ResolutionInfo *states = (ResolutionInfo *)calloc((size_t)num_states, sizeof(ResolutionInfo));

    compute_label_occurrences(knot, first_occ, second_occ);
    for (int r = 0; r < num_states; r++) compute_resolution_info(knot, first_occ, second_occ, r, &states[r]);

    int min_q = INT_MAX, max_q = INT_MIN;
    for (int r = 0; r < num_states; r++) {
        int h = pop_count_u32((uint32_t)r);
        int circles = states[r].n_circles;
        unsigned int basis_count = 1u << circles;
        for (unsigned int mask = 0; mask < basis_count; mask++) {
            int q = h + circles - 2 * pop_count_u32(mask) + n_plus - 2 * n_minus;
            if (q < min_q) min_q = q;
            if (q > max_q) max_q = q;
        }
    }

    int q_count = max_q - min_q + 1;
    BasisBucket **buckets = (BasisBucket **)malloc((size_t)(m + 1) * sizeof(BasisBucket *));
    int **ranks = (int **)malloc((size_t)m * sizeof(int *));
    int **betti = (int **)malloc((size_t)(m + 1) * sizeof(int *));
    
    int *torsion_counts = (int *)calloc((size_t)(m + 1) * q_count, sizeof(int));
    int *torsion_values = (int *)calloc((size_t)(m + 1) * q_count * 10, sizeof(int)); 

    for (int h = 0; h <= m; h++) {
        buckets[h] = (BasisBucket *)calloc((size_t)q_count, sizeof(BasisBucket));
        betti[h] = (int *)calloc((size_t)q_count, sizeof(int));
    }
    for (int h = 0; h < m; h++) ranks[h] = (int *)calloc((size_t)q_count, sizeof(int));

    for (int r = 0; r < num_states; r++) {
        int h = pop_count_u32((uint32_t)r);
        int circles = states[r].n_circles;
        unsigned int basis_count = 1u << circles;
        for (unsigned int mask = 0; mask < basis_count; mask++) {
            int q = h + circles - 2 * pop_count_u32(mask) + n_plus - 2 * n_minus;
            bucket_push(&buckets[h][q - min_q], r, mask);
        }
    }

    printf("left trefoil (n_+=%d, n_-=%d)\n", n_plus, n_minus);

    for (int h = 0; h < m; h++) {
        for (int q = min_q; q <= max_q; q++) {
            BasisBucket *source_bucket = &buckets[h][q - min_q];
            BasisBucket *target_bucket = &buckets[h + 1][q - min_q];
            Mat *differential = createMat(target_bucket->count, source_bucket->count);

            for (int col = 0; col < source_bucket->count; col++) {
                int src_resolution = source_bucket->resolutions[col];
                unsigned int src_mask = source_bucket->masks[col];
                for (int crossing = 0; crossing < m; crossing++) {
                    if ((src_resolution & (1 << crossing)) != 0) continue;
                    int dst_resolution = src_resolution | (1 << crossing);
                    int sign = jordan_wigner_sign(src_resolution, dst_resolution);
                    emit_edge_terms(&states[src_resolution], &states[dst_resolution], dst_resolution, src_mask, sign, target_bucket, differential, col);
                }
            }
            
            toSmithForm(differential);
            
            int current_rank = 0;
            int limit = differential->rows < differential->cols ? differential->rows : differential->cols;
            
            for (int i = 0; i < limit; i++) {
                int64_t val = differential->matrix[i][i];
                if (val != 0) {
                    current_rank++; 
                    if (val > 1 || val < -1) {
                        int64_t abs_val = val > 0 ? val : -val;
                        int idx = (h + 1) * q_count + (q - min_q);
                        int t_idx = torsion_counts[idx];
                        if (t_idx < 10) {
                            torsion_values[idx * 10 + t_idx] = (int)abs_val;
                            torsion_counts[idx]++;
                        }
                    }
                }
            }
            ranks[h][q - min_q] = current_rank;
            freeMat(differential);
        }
    }

    for (int h = 0; h <= m; h++) {
        for (int q = min_q; q <= max_q; q++) {
            int chain_rank = buckets[h][q - min_q].count;
            int rank_out = (h < m) ? ranks[h][q - min_q] : 0;
            int rank_in = (h > 0) ? ranks[h - 1][q - min_q] : 0;
            betti[h][q - min_q] = chain_rank - rank_out - rank_in;
        }
    }

    for (int h = 0; h <= m; h++) {
        for (int q = min_q; q <= max_q; q++) {
            int idx = h * q_count + (q - min_q);
            if (betti[h][q - min_q] != 0) {
                printf("rank H^{%d,%d} = %d\n", h - n_minus, q, betti[h][q - min_q]);
            }
            for (int i = 0; i < torsion_counts[idx]; i++) {
                printf("torsion H^{%d,%d} = Z_%d\n", h - n_minus, q, torsion_values[idx * 10 + i]);
            }
        }
    }
    
    print_bigraded_poincare(betti, m, min_q, max_q, n_minus, torsion_counts, torsion_values, q_count);

    for (int r = 0; r < num_states; r++) free_resolution_info(&states[r]);
    free(states); free(first_occ); free(second_occ);
    free(torsion_counts); free(torsion_values);
    for (int h = 0; h <= m; h++) {
        for (int q = 0; q < q_count; q++) bucket_free(&buckets[h][q]);
        free(buckets[h]); free(betti[h]);
    }
    for (int h = 0; h < m; h++) free(ranks[h]);
    free(buckets); free(ranks); free(betti); free_local_knot(knot);
    return 0;
}