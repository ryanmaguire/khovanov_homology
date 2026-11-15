#include "../KhMemory.h"
/*
 *  Komplex.c
 *
 *  Bar-Natan Khovanov Chain Complex with full Schur complement
 *  Gaussian elimination. Optional oracle scheduling is compiled only with
 *  KH_ENABLE_MORSE_ORACLE.
 *
 *  Purpose:
 *      Implements scanning-complex simplification.
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

/* ================================================================
 *  Construction / Destruction
 * ================================================================ */
Komplex *Komplex_create(int length) {
  if (length < 1) kh_fatal("invalid complex length");
  Komplex *k = (Komplex *)kh_malloc(sizeof(Komplex));
  k->length = length;
  if (length > 1) {
    k->differentials = (CobMatrix **)kh_calloc(length - 1, sizeof(CobMatrix *));
  } else {
    k->differentials = NULL;
  }
  k->chain_groups = (SmoothingColumn **)kh_calloc(length, sizeof(SmoothingColumn *));
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
    for (int i = 0; i < k->length; ++i)
      SmoothingColumn_free(k->chain_groups[i]);
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
  /* An LCCC is an isomorphism only if it has exactly one term whose
   * coefficient is a unit in the active coefficient system and whose
   * cobordism is itself an isomorphism. Over Z the only units are +/-1;
   * over F_p every nonzero coefficient is a unit. */
  if (lc->head != NULL && lc->head->next == NULL) {
    if (!LCCC_coeffIsUnit(lc->head->coeff)) return false;

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
  if (phi_inv == NULL)
    return false;
  /* ---- Step 3: Extract row B = D[target_row][*] (excluding source_col) ---- */
  /* And column C = D[*][source_col] (excluding target_row)              */
  /* Unpack the pivot row: B[j] = D[target_row][j] for j ≠ source_col */
  LCCC **B_row = CobMatrix_unpackRow(D, target_row);
  /* Unpack the pivot column: C[i] = D[i][source_col] for i ≠ target_row */
  int n_rows = D->target->n;
  int n_cols = D->source->n;
  LCCC **C_col = (LCCC **)kh_calloc((size_t)n_rows, sizeof(LCCC *));
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
       * Subtraction is implemented by coefficient-aware negation, so this is
       * valid both over Z and over the active prime field F_p:
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
  SmoothingColumn_free(k->chain_groups[chain_idx]);
  SmoothingColumn_free(k->chain_groups[chain_idx + 1]);
  k->chain_groups[chain_idx] = SmoothingColumn_retain(D->source);
  k->chain_groups[chain_idx + 1] = SmoothingColumn_retain(D->target);
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
      SmoothingColumn_free(d_prev->target);
      d_prev->target = SmoothingColumn_retain(k->chain_groups[chain_idx]);
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
      SmoothingColumn_free(d_next->source);
      d_next->source = SmoothingColumn_retain(k->chain_groups[chain_idx + 1]);
    }
  }
  return true;
}
/* ================================================================
 *  Oracle-Driven Reduction
 * ================================================================ */
#ifdef KH_ENABLE_MORSE_ORACLE
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
#endif
/* ================================================================
 *  Greedy Algebraic Reduction
 * ================================================================ */
int Komplex_greedyReduce(Komplex *k) {
  int reductions = 0;

  /*
   * Process each differential to a fixed point before advancing.
   *
   * Cancelling a pivot in d_i changes d_i by the Schur complement, while
   * d_{i-1} and d_{i+1} only lose a row/column. Those adjacent deletions do not
   * create new morphism values, so they cannot create a new isomorphism in a
   * differential that has already been exhausted. The only place a new pivot
   * can be created is the current d_i itself.
   *
   * This avoids restarting the entire chain complex at ci = 0 after every
   * successful cancellation.
   */
  for (int ci = 0; ci < k->length - 1; ci++) {
    while (true) {
      CobMatrix *D = k->differentials[ci];
      if (D == NULL)
        break;

      bool reduced_here = false;

      for (int tr = 0; tr < D->target->n && !reduced_here; tr++) {
        MatrixEntry *cur = D->entries[tr].head;
        while (cur != NULL) {
          int sc = cur->column_index;
          if (is_isomorphism(cur->value) && pivot_preserves_q(D, sc, tr) &&
              Komplex_blockReductionLemma(k, ci, sc, tr)) {
            reductions++;
            reduced_here = true;
            break; /* D changed; restart the search in this differential. */
          }
          cur = cur->next;
        }
      }

      if (!reduced_here)
        break;
    }
  }

  return reductions;
}
/*
 * A generator whose cap has c circles splits into 2^c generators. Circles
 * are removed last-first; variant j has bit (c-1-s) set when step s took the
 * q-1 summand, so its grading is base + c - 2*popcount(j).
 *
 * expand_entry replaces one matrix entry by its 2^c delooped entries.
 * caps[s] is the generator's cap with s circles removed. incoming selects
 * the differential into the column (cups on the bottom) or out of it (caps
 * on the top). The input value is borrowed; outputs are owned by the caller.
 */
static LCCC **expand_entry(LCCC *value, Cap **caps, int cycles, bool incoming) {
  size_t count = (size_t)1 << cycles;
  LCCC **cur = (LCCC **)kh_malloc(count * sizeof(LCCC *));
  LCCC **next = (LCCC **)kh_malloc(count * sizeof(LCCC *));
  size_t have = 1;
  cur[0] = value;
  for (int s = 1; s <= cycles; s++) {
    for (size_t v = 0; v < have; v++) {
      if (incoming) {
        next[2 * v] = LCCC_cupOnBottom(cur[v], caps[s], true);
        next[2 * v + 1] = LCCC_cupOnBottom(cur[v], caps[s], false);
      } else {
        next[2 * v] = LCCC_capOffTop(cur[v], caps[s], false);
        next[2 * v + 1] = LCCC_capOffTop(cur[v], caps[s], true);
      }
      if (s > 1) LCCC_free(cur[v]);
    }
    LCCC **tmp = cur;
    cur = next;
    next = tmp;
    have *= 2;
  }
  free(next);
  return cur;
}

/* Entries of generators without circles are moved, not copied, so callers
 * must have reduced the differentials beforehand. Expanded entries come out
 * reduced. */
void Komplex_deloop(Komplex *k) {
  for (int h = 0; h < k->length; h++) {
    SmoothingColumn *col = k->chain_groups[h];
    if (col == NULL || col->n == 0)
      continue;

    int old_n = col->n;
    bool any = false;
    int64_t total = 0;
    int max_cycles = 0;
    for (int i = 0; i < old_n; i++) {
      if (col->smoothings[i] == NULL)
        kh_fatal("missing cap in delooping");
      int c = col->smoothings[i]->ncycles;
      if (c > 30) kh_fatal("too many circles to deloop");
      if (c > 0) any = true;
      if (c > max_cycles) max_cycles = c;
      total += (int64_t)1 << c;
    }
    if (!any)
      continue;

    int new_n = kh_int(total);
    int *start = (int *)kh_malloc((size_t)old_n * sizeof(int));
    /* caps[i][s]: cap of generator i with s circles removed (borrowed at s=0). */
    Cap ***caps = (Cap ***)kh_calloc((size_t)old_n, sizeof(Cap **));

    SmoothingColumn *new_col = SmoothingColumn_create();
    new_col->n = new_n;
    new_col->numbers = (int *)kh_malloc((size_t)new_n * sizeof(int));
    new_col->smoothings = (Cap **)kh_malloc((size_t)new_n * sizeof(Cap *));

    int at = 0;
    for (int i = 0; i < old_n; i++) {
      Cap *base = col->smoothings[i];
      int c = base->ncycles;
      start[i] = at;
      if (c == 0) {
        new_col->numbers[at] = col->numbers[i];
        new_col->smoothings[at] = Cap_retain(base);
        at++;
        continue;
      }
      caps[i] = (Cap **)kh_malloc((size_t)(c + 1) * sizeof(Cap *));
      caps[i][0] = base;
      for (int s = 1; s <= c; s++) {
        caps[i][s] = Cap_removeCycle(caps[i][s - 1]);
        if (caps[i][s] == NULL)
          kh_fatal("failed to remove circle in delooping");
      }
      int variants = 1 << c;
      for (int j = 0; j < variants; j++) {
        int minus = 0;
        for (int b = j; b; b >>= 1) minus += b & 1;
        new_col->numbers[at] = kh_int((int64_t)col->numbers[i] + c - 2 * minus);
        new_col->smoothings[at] = Cap_retain(caps[i][c]);
        at++;
      }
    }

    if (h > 0 && k->differentials[h - 1] != NULL) {
      CobMatrix *old_in = k->differentials[h - 1];
      CobMatrix *new_in = CobMatrix_create(old_in->source, new_col, true);
      for (int r = 0; r < old_n; r++) {
        int c = col->smoothings[r]->ncycles;
        for (MatrixEntry *e = old_in->entries[r].head; e != NULL; e = e->next) {
          if (c == 0) {
            CobMatrix_putEntry(new_in, start[r], e->column_index, e->value);
            e->value = NULL;
            continue;
          }
          LCCC **out = expand_entry(e->value, caps[r], c, true);
          int variants = 1 << c;
          for (int j = 0; j < variants; j++)
            CobMatrix_putEntry(new_in, start[r] + j, e->column_index, out[j]);
          free(out);
        }
      }
      k->differentials[h - 1] = new_in;
      CobMatrix_free(old_in);
    }

    if (h < k->length - 1 && k->differentials[h] != NULL) {
      CobMatrix *old_out = k->differentials[h];
      CobMatrix *new_out = CobMatrix_create(new_col, old_out->target, true);
      for (int r = 0; r < old_out->target->n; r++) {
        for (MatrixEntry *e = old_out->entries[r].head; e != NULL; e = e->next) {
          int src = e->column_index;
          int c = col->smoothings[src]->ncycles;
          if (c == 0) {
            CobMatrix_putEntry(new_out, r, start[src], e->value);
            e->value = NULL;
            continue;
          }
          LCCC **out = expand_entry(e->value, caps[src], c, false);
          int variants = 1 << c;
          for (int j = 0; j < variants; j++)
            CobMatrix_putEntry(new_out, r, start[src] + j, out[j]);
          free(out);
        }
      }
      k->differentials[h] = new_out;
      CobMatrix_free(old_out);
    }

    for (int i = 0; i < old_n; i++) {
      if (caps[i] == NULL) continue;
      for (int s = 1; s <= col->smoothings[i]->ncycles; s++)
        Cap_free(caps[i][s]);
      free(caps[i]);
    }
    free(caps);
    free(start);

    k->chain_groups[h] = new_col;
    SmoothingColumn_free(col);
  }
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
      return false;
    CobMatrix_reduce(composed);
    if (!CobMatrix_isZero(composed)) {
      CobMatrix_free(composed);
      return false;
    }
    CobMatrix_free(composed);
  }
  return true;
}
