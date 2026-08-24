/*
 *  Komplex.c
 *
 *  Bar-Natan Khovanov Chain Complex with full Schur complement
 *  Gaussian elimination, guided by the topological Morse oracle.
 *
 *  Purpose:
 *      Implements the algebraic Phase 2 of the FPT algorithm.
 *      The blockReductionLemma performs the complete Schur complement
 *      update D' = D - C·φ⁻¹·B on the differential, and propagates
 *      generator removal to adjacent differentials.
 *
 *  Compile:
 *      gcc -Wall -Wextra -g -c Komplex.c -I.
 */
#include "Komplex.h"
#include "LCCC.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static void CobMatrix_free_with_known_row_count(CobMatrix *m, int row_count) {
  if (m == NULL)
    return;

  if (m->entries != NULL) {
    for (int i = 0; i < row_count; i++) {
      MatrixEntry *cur = m->entries[i].head;
      while (cur != NULL) {
        MatrixEntry *next = cur->next;

        LCCC_free(cur->value);
        free(cur);

        cur = next;
      }
    }

    free(m->entries);
  }

  free(m);
}

/* ================================================================
 *  Construction / Destruction
 * ================================================================ */
Komplex *Komplex_create(int length) {
  Komplex *k = (Komplex *)malloc(sizeof(Komplex));
  k->length = length;
  if (length > 1) {
    k->differentials = (CobMatrix **)calloc(length - 1, sizeof(CobMatrix *));
  } else {
    k->differentials = NULL;
  }
  k->chain_groups = (SmoothingColumn **)calloc(length, sizeof(SmoothingColumn *));
  return k;
}
void Komplex_free(Komplex *k) {
  if (!k)
    return;
  for (int i = 0; i < k->length - 1; i++) {
    if (k->differentials[i]) {
      CobMatrix_free(k->differentials[i]);
    }
  }
  free(k->differentials);
  
  if (k->chain_groups) {
    /* Ownership of SmoothingColumns is tricky; if they are managed by chain_groups, free them.
     * But they might be aliased by CobMatrix source/target.
     * For now we just free the array. */
    free(k->chain_groups);
  }
  free(k);
}
/* ================================================================
 *  Block Reduction Lemma — Full Schur Complement
 *
 *  Given: D = differentials[chain_idx] with pivot φ = D[target_row][source_col]
 *
 *  We write D in block form:
 *
 *              source_col    rest
 *  target_row [   φ    |    B   ]
 *  rest       [   C    |    D₀  ]
 *
 *  The reduced differential is:
 *      D' = D₀ - C · φ⁻¹ · B
 *
 *  Additionally:
 *  - d_{chain_idx-1} loses its row corresponding to source_col
 *    (since source_col is a target index in d_{chain_idx-1})
 *  - d_{chain_idx+1} loses its column corresponding to target_row
 *    (since target_row is a source index in d_{chain_idx+1})
 * ================================================================ */
static bool is_isomorphism(LCCC *lc) {
  if (lc == NULL || LCCC_isZero(lc)) return false;
  /* An LCCC is an isomorphism only if it has exactly one term
     with coefficient ±1 and that term is an isomorphism. */
  if (lc->head != NULL && lc->head->next == NULL) {
    if (lc->head->coeff != 1 && lc->head->coeff != -1) return false;

    CannedCobordism *cc = lc->head->cobordism;
    if (cc != NULL && cc->isIsomorphism != NULL) {
      return cc->isIsomorphism(cc);
    }
  }
  return false;
}
static bool pivot_preserves_q(CobMatrix *D, int source_col, int target_row) {
  if (D == NULL || D->source == NULL || D->target == NULL)
    return false;
  if (D->source->numbers == NULL || D->target->numbers == NULL)
    return false;
  if (source_col < 0 || source_col >= D->source->n)
    return false;
  if (target_row < 0 || target_row >= D->target->n)
    return false;

  return D->source->numbers[source_col] == D->target->numbers[target_row];
}
bool Komplex_blockReductionLemma(Komplex *k, int chain_idx, int source_col,
                                 int target_row) {
  if (chain_idx < 0 || chain_idx >= k->length - 1)
    return false;
  CobMatrix *D = k->differentials[chain_idx];
  if (!D)
    return false;
  /* Validate indices */
  if (source_col < 0 || source_col >= D->source->n)
    return false;
  if (target_row < 0 || target_row >= D->target->n)
    return false;
  /* ---- Step 1: Extract the pivot φ = D[target_row][source_col] ---- */
  LCCC *phi = NULL;
  {
    MatrixEntry *cur = D->entries[target_row].head;
    while (cur != NULL) {
      if (cur->column_index == source_col) {
        phi = cur->value;
        break;
      }
      cur = cur->next;
    }
  }
  /* Verify that the selected pivot is an invertible term. */
  if (!is_isomorphism(phi))
    return false;
  if (!pivot_preserves_q(D, source_col, target_row))
    return false;
  /* ---- Step 2: Compute φ⁻¹ ---- */
  LCCC *phi_inv = LCCC_invert(phi);
  /* ---- Step 3: Extract row B = D[target_row][*] (excluding source_col) ---- */
  /* And column C = D[*][source_col] (excluding target_row)              */
  /* Unpack the pivot row: B[j] = D[target_row][j] for j ≠ source_col */
  LCCC **B_row = CobMatrix_unpackRow(D, target_row);
  /* Unpack the pivot column: C[i] = D[i][source_col] for i ≠ target_row */
  int n_rows = D->target->n;
  int n_cols = D->source->n;
  LCCC **C_col = (LCCC **)calloc((size_t)n_rows, sizeof(LCCC *));
  for (int i = 0; i < n_rows; i++) {
    MatrixEntry *cur = D->entries[i].head;
    while (cur != NULL) {
      if (cur->column_index == source_col) {
        C_col[i] = cur->value;
        break;
      }
      cur = cur->next;
    }
  }
  /* ---- Step 4: Schur complement: D[i][j] -= C[i] · φ⁻¹ · B[j] ---- */
  for (int i = 0; i < n_rows; i++) {
    if (i == target_row) continue;
    if (C_col[i] == NULL || LCCC_isZero(C_col[i])) continue;
    /* Compute C[i] · φ⁻¹ */
    LCCC *ci_phi_inv = LCCC_compose(C_col[i], phi_inv);
    if (ci_phi_inv == NULL || LCCC_isZero(ci_phi_inv)) {
      LCCC_free(ci_phi_inv);
      continue;
    }
    for (int j = 0; j < n_cols; j++) {
      if (j == source_col) continue;
      if (B_row[j] == NULL || LCCC_isZero(B_row[j])) continue;
      /* Compute ci_phi_inv · B[j] = C[i] · φ⁻¹ · B[j] */
      LCCC *correction = LCCC_compose(ci_phi_inv, B_row[j]);
      if (correction == NULL || LCCC_isZero(correction)) {
        LCCC_free(correction);
        continue;
      }
      /*
      * Over Z, subtraction is implemented by adding the negated correction:
      * D[i][j] -= correction.
      */
      LCCC *neg_correction = LCCC_negate(correction);
      CobMatrix_addEntry(D, i, j, neg_correction);
      LCCC_free(correction);
      /* Note: neg_correction ownership is transferred to addEntry */
    }
    LCCC_free(ci_phi_inv);
  }
  free(B_row);
  free(C_col);
  LCCC_free(phi_inv);
  /* ---- Step 5: Remove source_col and target_row from D ---- */
  CobMatrix *removed_col = CobMatrix_extractColumn(D, source_col);
  CobMatrix *removed_row = CobMatrix_extractRow(D, target_row);
  CobMatrix_free(removed_col);
  CobMatrix_free(removed_row);
  k->chain_groups[chain_idx] = D->source;
  k->chain_groups[chain_idx + 1] = D->target;
  /* ---- Step 6: Update adjacent differential d_{chain_idx - 1} ---- */
  /*
   * d_{prev} maps C_{chain_idx-1} → C_{chain_idx}.
   * d_{prev}->target = C_{chain_idx} generators.
   * We removed source_col from C_{chain_idx}, so remove that row
   * from d_{prev}.
   */
  if (chain_idx > 0 && k->differentials[chain_idx - 1] != NULL) {
    CobMatrix *d_prev = k->differentials[chain_idx - 1];
    if (source_col < d_prev->target->n) {
      CobMatrix *prev_removed = CobMatrix_extractRow(d_prev, source_col);
      CobMatrix_free(prev_removed);
      d_prev->target = k->chain_groups[chain_idx];
    }
  }
  /* ---- Step 7: Update adjacent differential d_{chain_idx + 1} ---- */
  /*
   * d_{next} maps C_{chain_idx+1} → C_{chain_idx+2}.
   * d_{next}->source = C_{chain_idx+1} generators.
   * We removed target_row from C_{chain_idx+1}, so remove that
   * column from d_{next}.
   */
  if (chain_idx + 1 < k->length - 1 && k->differentials[chain_idx + 1] != NULL) {
    CobMatrix *d_next = k->differentials[chain_idx + 1];
    if (target_row < d_next->source->n) {
      CobMatrix *next_removed = CobMatrix_extractColumn(d_next, target_row);
      CobMatrix_free(next_removed);
      d_next->source = k->chain_groups[chain_idx + 1];
    }
  }
  return true;
}
/* ================================================================
 *  Oracle-Driven Reduction
 * ================================================================ */
int Komplex_reduce_with_oracle(Komplex *k, const CollapseSchedule *schedule) {
  int reductions = 0;
  for (int step = 0; step < schedule->count; step++) {
    int ci = schedule->chain_idx[step];
    int sc = schedule->source_col[step];
    int tr = schedule->target_row[step];
    /* Skip invalid entries (unmapped pairs) */
    if (ci < 0 || sc < 0 || tr < 0)
      continue;
    if (Komplex_blockReductionLemma(k, ci, sc, tr)) {
      reductions++;
    }
  }
  return reductions;
}
/* ================================================================
 *  Greedy Algebraic Reduction
 * ================================================================ */
int Komplex_greedyReduce(Komplex *k) {
  int reductions = 0;
  bool reduced;
  do {
    reduced = false;
    for (int ci = 0; ci < k->length - 1; ci++) {
      CobMatrix *D = k->differentials[ci];
      if (!D) continue;
      
      /* Search for an isomorphism in D */
      for (int tr = 0; tr < D->target->n; tr++) {
        MatrixEntry *cur = D->entries[tr].head;
        while (cur != NULL) {
          int sc = cur->column_index;
          if (is_isomorphism(cur->value) && pivot_preserves_q(D, sc, tr)) {
            if (Komplex_blockReductionLemma(k, ci, sc, tr)) {
              reductions++;
              reduced = true;
              break; /* Break out of while loop, D has changed */
            }
          }
          cur = cur->next;
        }
        if (reduced) break; /* Break out of target loop */
      }
      if (reduced) break; /* Restart the full scan from ci=0 because adjacent matrices changed */
    }
  } while (reduced);
  
  return reductions;
}
void Komplex_deloop(Komplex *k) {
  bool reduced;
  do {
    reduced = false;

    for (int h = 0; h < k->length; h++) {
      SmoothingColumn *col = k->chain_groups[h];
      if (col == NULL)
        continue;

      int split_idx = -1;
      for (int i = 0; i < col->n; i++) {
        if (col->smoothings[i] != NULL && col->smoothings[i]->ncycles > 0) {
          split_idx = i;
          break;
        }
      }
      if (split_idx == -1)
        continue;

      int old_n = col->n;
      int new_n = old_n + 1;
      int *new_numbers = (int *)malloc((size_t)new_n * sizeof(int));
      Cap **new_smoothings = (Cap **)malloc((size_t)new_n * sizeof(Cap *));
      if (new_numbers == NULL || new_smoothings == NULL) {
        free(new_numbers);
        free(new_smoothings);
        return;
      }

      memcpy(new_numbers, col->numbers, (size_t)old_n * sizeof(int));
      memcpy(new_smoothings, col->smoothings, (size_t)old_n * sizeof(Cap *));
      free(col->numbers);
      free(col->smoothings);

      col->n = new_n;
      col->numbers = new_numbers;
      col->smoothings = new_smoothings;

      Cap *base_cap = col->smoothings[split_idx];
      if (base_cap == NULL)
        return;
      Cap *v_plus_cap = Cap_removeCycle(base_cap);
      Cap *v_minus_cap = Cap_removeCycle(base_cap);
      if (v_plus_cap == NULL || v_minus_cap == NULL)
        return;

      col->smoothings[split_idx] = v_plus_cap;
      col->smoothings[old_n] = v_minus_cap;

      int base_shift = col->numbers[split_idx];
      col->numbers[split_idx] = base_shift + 1;
      col->numbers[old_n] = base_shift - 1;

      if (h > 0 && k->differentials[h - 1] != NULL) {
        CobMatrix *old_d_in = k->differentials[h - 1];
        CobMatrix *new_d_in = CobMatrix_create(old_d_in->source, col, true);

        for (int r = 0; r < old_n; r++) {
          MatrixEntry *cur = old_d_in->entries[r].head;

          while (cur != NULL) {
            int c = cur->column_index;

            if (r == split_idx) {
              LCCC *lccc_plus =
                  LCCC_cupOnBottom(cur->value, v_plus_cap, true);
              LCCC *lccc_minus =
                  LCCC_cupOnBottom(cur->value, v_minus_cap, false);

              CobMatrix_addEntry(new_d_in, split_idx, c, lccc_plus);
              CobMatrix_addEntry(new_d_in, old_n, c, lccc_minus);
            } else {
              CobMatrix_addEntry(new_d_in, r, c, LCCC_clone(cur->value));
            }

            cur = cur->next;
          }
        }

        CobMatrix_reduce(new_d_in);

        k->differentials[h - 1] = new_d_in;

        CobMatrix_free_with_known_row_count(old_d_in, old_n);
      }

      if (h < k->length - 1 && k->differentials[h] != NULL) {
        CobMatrix *old_d_out = k->differentials[h];
        CobMatrix *new_d_out = CobMatrix_create(col, old_d_out->target, true);

        for (int r = 0; r < old_d_out->target->n; r++) {
          MatrixEntry *cur = old_d_out->entries[r].head;

          while (cur != NULL) {
            int c = cur->column_index;

            if (c == split_idx) {
              LCCC *lccc_plus =
                  LCCC_capOffTop(cur->value, v_plus_cap, false);
              LCCC *lccc_minus =
                  LCCC_capOffTop(cur->value, v_minus_cap, true);

              CobMatrix_addEntry(new_d_out, r, split_idx, lccc_plus);
              CobMatrix_addEntry(new_d_out, r, old_n, lccc_minus);
            } else {
              CobMatrix_addEntry(new_d_out, r, c, LCCC_clone(cur->value));
            }

            cur = cur->next;
          }
        }

        CobMatrix_reduce(new_d_out);

        k->differentials[h] = new_d_out;
        CobMatrix_free(old_d_out);
      }
      for (int di = 0; di < k->length - 1; di++) {
        if (k->differentials[di] != NULL) {
          CobMatrix_reduce(k->differentials[di]);
        }
      }

      reduced = true;
      break;
    }
  } while (reduced);
}
/* ================================================================
 *  Chain Complex Verification: d² = 0
 * ================================================================ */
bool Komplex_verify_d_squared(const Komplex *k) {
  if (k->length <= 2)
    return true; /* Nothing to check with 0 or 1 differentials */
  for (int i = 0; i < k->length - 2; i++) {
    CobMatrix *d_i = k->differentials[i];
    CobMatrix *d_ip1 = k->differentials[i + 1];
    if (d_i == NULL || d_ip1 == NULL)
      continue;
    /* Compute d_{i+1} ∘ d_i */
    CobMatrix *composed = CobMatrix_compose(d_ip1, d_i);
    if (composed == NULL)
      continue;
    CobMatrix_reduce(composed);
    if (!CobMatrix_isZero(composed)) {
      CobMatrix_free(composed);
      return false;
    }
    CobMatrix_free(composed);
  }
  return true;
}
