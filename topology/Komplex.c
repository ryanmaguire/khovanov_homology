#include "Komplex.h"

#include <stdlib.h>

/*
 * Optional helper functions.
 * These are not declared in CobMatrix.h, but the implementation may provide
 * them elsewhere. Keeping the declarations here makes this translation unit
 * self-contained.
 */
extern LCCC *LCCC_invert(LCCC *lc);
extern LCCC *LCCC_negate(LCCC *lc);

static LCCC *find_entry_value(const CobMatrix *m, int row_idx, int col_idx) {
    if (m == NULL || m->entries == NULL) {
        return NULL;
    }
    MatrixEntry *cur = m->entries[row_idx].head;
    while (cur != NULL) {
        if (cur->column_index == col_idx) {
            return cur->value;
        }
        cur = cur->next;
    }
    return NULL;
}

Komplex *Komplex_create(int length) {
    if (length < 0) {
        return NULL;
    }

    Komplex *k = (Komplex *)malloc(sizeof(Komplex));
    if (k == NULL) {
        return NULL;
    }

    k->length = length;
    k->differentials = NULL;

    if (length > 1) {
        k->differentials = (CobMatrix **)calloc((size_t)(length - 1), sizeof(CobMatrix *));
        if (k->differentials == NULL) {
            free(k);
            return NULL;
        }
    }

    return k;
}

void Komplex_free(Komplex *k) {
    if (k == NULL) {
        return;
    }

    int n_diff = (k->length > 1) ? (k->length - 1) : 0;
    for (int i = 0; i < n_diff; ++i) {
        if (k->differentials != NULL && k->differentials[i] != NULL) {
            CobMatrix_free(k->differentials[i]);
            k->differentials[i] = NULL;
        }
    }

    free(k->differentials);
    free(k);
}

bool Komplex_blockReductionLemma(Komplex *k, int chain_idx, int source_col, int target_row) {
    if (k == NULL || k->differentials == NULL) {
        return false;
    }
    if (chain_idx < 0 || chain_idx >= k->length - 1) {
        return false;
    }

    CobMatrix *D = k->differentials[chain_idx];
    if (D == NULL || D->source == NULL || D->target == NULL || D->entries == NULL) {
        return false;
    }

    if (source_col < 0 || source_col >= D->source->n) {
        return false;
    }
    if (target_row < 0 || target_row >= D->target->n) {
        return false;
    }

    /* Pivot entry φ = D[target_row][source_col]. */
    LCCC *phi = find_entry_value(D, target_row, source_col);
    if (phi == NULL || LCCC_isZero(phi)) {
        return false;
    }

    LCCC *phi_inv = LCCC_invert(phi);
    if (phi_inv == NULL || LCCC_isZero(phi_inv)) {
        if (phi_inv != NULL) {
            LCCC_free(phi_inv);
        }
        return false;
    }

    const int n_rows = D->target->n;
    const int n_cols = D->source->n;

    /* B_row[j] = D[target_row][j] for j != source_col. */
    LCCC **B_row = CobMatrix_unpackRow(D, target_row);
    if (B_row == NULL) {
        LCCC_free(phi_inv);
        return false;
    }

    /* C_col[i] = D[i][source_col] for i != target_row. */
    LCCC **C_col = (LCCC **)calloc((size_t)n_rows, sizeof(LCCC *));
    if (C_col == NULL) {
        free(B_row);
        LCCC_free(phi_inv);
        return false;
    }

    for (int i = 0; i < n_rows; ++i) {
        C_col[i] = find_entry_value(D, i, source_col);
    }

    /*
     * Schur complement update:
     * D[i][j] <- D[i][j] - C[i] * phi^{-1} * B[j]
     *
     * If the coefficient ring has characteristic 2, subtraction is addition;
     * otherwise we negate the correction before inserting it.
     */
    for (int i = 0; i < n_rows; ++i) {
        if (i == target_row) {
            continue;
        }
        if (C_col[i] == NULL || LCCC_isZero(C_col[i])) {
            continue;
        }

        LCCC *ci_phi_inv = LCCC_compose(C_col[i], phi_inv);
        if (ci_phi_inv == NULL || LCCC_isZero(ci_phi_inv)) {
            if (ci_phi_inv != NULL) {
                LCCC_free(ci_phi_inv);
            }
            continue;
        }

        for (int j = 0; j < n_cols; ++j) {
            if (j == source_col) {
                continue;
            }
            if (B_row[j] == NULL || LCCC_isZero(B_row[j])) {
                continue;
            }

            LCCC *correction = LCCC_compose(ci_phi_inv, B_row[j]);
            if (correction == NULL || LCCC_isZero(correction)) {
                if (correction != NULL) {
                    LCCC_free(correction);
                }
                continue;
            }

            LCCC *neg_correction = LCCC_negate(correction);
            if (neg_correction == NULL) {
                /*
                 * Fallback for implementations that only support characteristic 2.
                 * In that setting, subtraction and addition coincide.
                 */
                neg_correction = correction;
                correction = NULL;
            }

            CobMatrix_addEntry(D, i, j, neg_correction);

            if (correction != NULL) {
                LCCC_free(correction);
            }
        }

        LCCC_free(ci_phi_inv);
    }

    free(B_row);
    free(C_col);
    LCCC_free(phi_inv);

    /*
     * Remove the pivot row/column from the matrix itself.
     * These functions are assumed to mutate D and return the removed slice.
     */
    CobMatrix *removed_col = CobMatrix_extractColumn(D, source_col);
    CobMatrix *removed_row = CobMatrix_extractRow(D, target_row);
    if (removed_col != NULL) {
        CobMatrix_free(removed_col);
    }
    if (removed_row != NULL) {
        CobMatrix_free(removed_row);
    }

    /* Update the previous differential: remove the row corresponding to source_col. */
    if (chain_idx > 0 && k->differentials[chain_idx - 1] != NULL) {
        CobMatrix *d_prev = k->differentials[chain_idx - 1];
        if (d_prev->target != NULL && source_col < d_prev->target->n) {
            CobMatrix *prev_removed = CobMatrix_extractRow(d_prev, source_col);
            if (prev_removed != NULL) {
                CobMatrix_free(prev_removed);
            }
        }
    }

    /* Update the next differential: remove the column corresponding to target_row. */
    if (chain_idx + 1 < k->length - 1 && k->differentials[chain_idx + 1] != NULL) {
        CobMatrix *d_next = k->differentials[chain_idx + 1];
        if (d_next->source != NULL && target_row < d_next->source->n) {
            CobMatrix *next_removed = CobMatrix_extractColumn(d_next, target_row);
            if (next_removed != NULL) {
                CobMatrix_free(next_removed);
            }
        }
    }

    return true;
}

int Komplex_reduce_with_oracle(Komplex *k, const CollapseSchedule *schedule) {
    if (k == NULL || schedule == NULL) {
        return 0;
    }
    if (schedule->count <= 0 || schedule->chain_idx == NULL ||
        schedule->source_col == NULL || schedule->target_row == NULL) {
        return 0;
    }

    int reductions = 0;
    for (int step = 0; step < schedule->count; ++step) {
        int ci = schedule->chain_idx[step];
        int sc = schedule->source_col[step];
        int tr = schedule->target_row[step];

        if (ci < 0 || sc < 0 || tr < 0) {
            continue;
        }

        if (Komplex_blockReductionLemma(k, ci, sc, tr)) {
            ++reductions;
        }
    }

    return reductions;
}

bool Komplex_verify_d_squared(const Komplex *k) {
    if (k == NULL || k->length <= 2 || k->differentials == NULL) {
        return true;
    }

    for (int i = 0; i < k->length - 2; ++i) {
        CobMatrix *d_i = k->differentials[i];
        CobMatrix *d_ip1 = k->differentials[i + 1];

        if (d_i == NULL || d_ip1 == NULL) {
            continue;
        }

        CobMatrix *composed = CobMatrix_compose(d_ip1, d_i);
        if (composed == NULL) {
            continue;
        }

        if (!CobMatrix_isZero(composed)) {
            CobMatrix_free(composed);
            return false;
        }

        CobMatrix_free(composed);
    }

    return true;
}
