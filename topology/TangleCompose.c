#include "TangleCompose.h"
#include "CannedCobordismImpl.h"
#include "LCCC.h"
#include <stdlib.h>

static SmoothingColumn *tensor_column(SmoothingColumn *a, int astart,
                                      SmoothingColumn *b, int bstart,
                                      int join_count) {
  if (a == NULL || b == NULL) return NULL;

  SmoothingColumn *out = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  if (out == NULL) return NULL;

  out->n = a->n * b->n;
  out->numbers = out->n > 0 ? (int *)malloc((size_t)out->n * sizeof(int)) : NULL;
  out->smoothings =
      out->n > 0 ? (Cap **)malloc((size_t)out->n * sizeof(Cap *)) : NULL;
  if (out->n > 0 && (out->numbers == NULL || out->smoothings == NULL)) {
    free(out->numbers);
    free(out->smoothings);
    free(out);
    return NULL;
  }

  int index = 0;
  for (int i = 0; i < a->n; i++) {
    for (int j = 0; j < b->n; j++) {
      out->numbers[index] = a->numbers[i] + b->numbers[j];
      out->smoothings[index] =
          Cap_compose(a->smoothings[i], astart, b->smoothings[j], bstart,
                      join_count, NULL);
      if (out->smoothings[index] == NULL) {
        for (int k = 0; k < index; k++) Cap_free(out->smoothings[k]);
        free(out->numbers);
        free(out->smoothings);
        free(out);
        return NULL;
      }
      index++;
    }
  }
  return out;
}

static LCCC *compose_lccc_with_identity(const LCCC *lc, int lc_start,
                                        Cap *identity_cap, int id_start,
                                        int join_count, int coefficient_sign,
                                        bool lc_on_left) {
  if (lc == NULL) return LCCC_createZero();

  LCCC *result = LCCC_createZero();
  if (result == NULL) return NULL;

  CannedCobordism *identity = CannedCobordismImpl_isomorphism(identity_cap);
  if (identity == NULL) {
    LCCC_free(result);
    return NULL;
  }

  for (LCCCTerm *term = lc->head; term != NULL; term = term->next) {
    CannedCobordism *composed = NULL;
    if (lc_on_left) {
      composed = CannedCobordism_compose_partial(
          term->cobordism, lc_start, identity, id_start, join_count);
    } else {
      composed = CannedCobordism_compose_partial(
          identity, id_start, term->cobordism, lc_start, join_count);
    }
    if (composed == NULL) {
      CannedCobordism_free(identity);
      LCCC_free(result);
      return NULL;
    }

    LCCC *single =
        LCCC_createSingle(composed, coefficient_sign * term->coeff);
    if (single == NULL) {
      CannedCobordism_free(identity);
      LCCC_free(result);
      return NULL;
    }
    LCCC *sum = LCCC_add(result, single);
    LCCC_free(result);
    LCCC_free(single);
    if (sum == NULL) {
      CannedCobordism_free(identity);
      return NULL;
    }
    result = sum;
  }

  CannedCobordism_free(identity);
  return result;
}

Komplex *Komplex_compose_partial_tangles(Komplex *left, int left_start,
                                         Komplex *right, int right_start,
                                         int join_count) {
  if (left == NULL || right == NULL || join_count < 0) return NULL;

  int new_length = left->length + right->length - 1;
  Komplex *out = Komplex_create(new_length);
  if (out == NULL) return NULL;

  for (int h = 0; h < new_length; h++) {
    int total = 0;
    for (int i = 0; i < left->length; i++) {
      int j = h - i;
      if (j >= 0 && j < right->length)
        total += left->chain_groups[i]->n * right->chain_groups[j]->n;
    }

    SmoothingColumn *column = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
    if (column == NULL) return NULL;
    column->n = total;
    column->numbers = total > 0 ? (int *)malloc((size_t)total * sizeof(int)) : NULL;
    column->smoothings =
        total > 0 ? (Cap **)malloc((size_t)total * sizeof(Cap *)) : NULL;
    if (total > 0 && (column->numbers == NULL || column->smoothings == NULL))
      return NULL;

    int index = 0;
    for (int i = 0; i < left->length; i++) {
      int j = h - i;
      if (j < 0 || j >= right->length) continue;
      SmoothingColumn *piece =
          tensor_column(left->chain_groups[i], left_start,
                        right->chain_groups[j], right_start, join_count);
      if (piece == NULL) return NULL;
      for (int k = 0; k < piece->n; k++) {
        column->numbers[index] = piece->numbers[k];
        column->smoothings[index] = piece->smoothings[k];
        index++;
      }
      free(piece->numbers);
      free(piece->smoothings);
      free(piece);
    }
    out->chain_groups[h] = column;
  }

  for (int h = 0; h < new_length - 1; h++) {
    CobMatrix *d =
        CobMatrix_create(out->chain_groups[h], out->chain_groups[h + 1], true);
    if (d == NULL) return NULL;

    int row_offset = 0;
    for (int i_out = 0; i_out < left->length; i_out++) {
      int j_out = h + 1 - i_out;
      if (j_out < 0 || j_out >= right->length) continue;

      int col_offset = 0;
      for (int i_in = 0; i_in < left->length; i_in++) {
        int j_in = h - i_in;
        if (j_in < 0 || j_in >= right->length) continue;

        if (i_out == i_in + 1 && j_out == j_in) {
          CobMatrix *dl = left->differentials[i_in];
          for (int r = 0; r < dl->target->n; r++) {
            for (MatrixEntry *entry = dl->entries[r].head; entry != NULL;
                 entry = entry->next) {
              for (int b = 0; b < right->chain_groups[j_in]->n; b++) {
                int row = row_offset + r * right->chain_groups[j_in]->n + b;
                int col = col_offset +
                          entry->column_index * right->chain_groups[j_in]->n + b;
                LCCC *value = compose_lccc_with_identity(
                    entry->value, left_start,
                    right->chain_groups[j_in]->smoothings[b], right_start,
                    join_count, 1, true);
                if (value == NULL) return NULL;
                CobMatrix_addEntry(d, row, col, value);
              }
            }
          }
        }

        if (i_out == i_in && j_out == j_in + 1) {
          CobMatrix *dr = right->differentials[j_in];
          int sign = (i_in & 1) ? -1 : 1;
          for (int r = 0; r < dr->target->n; r++) {
            for (MatrixEntry *entry = dr->entries[r].head; entry != NULL;
                 entry = entry->next) {
              for (int a = 0; a < left->chain_groups[i_in]->n; a++) {
                int row = row_offset + a * dr->target->n + r;
                int col = col_offset + a * dr->source->n + entry->column_index;
                LCCC *value = compose_lccc_with_identity(
                    entry->value, right_start,
                    left->chain_groups[i_in]->smoothings[a], left_start,
                    join_count, sign, false);
                if (value == NULL) return NULL;
                CobMatrix_addEntry(d, row, col, value);
              }
            }
          }
        }

        col_offset += left->chain_groups[i_in]->n * right->chain_groups[j_in]->n;
      }

      row_offset += left->chain_groups[i_out]->n * right->chain_groups[j_out]->n;
    }
    out->differentials[h] = d;
  }

  return out;
}
