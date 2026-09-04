#include "TangleKomplex.h"
#include "Komplex.h"
#include "Cap.h"
#include "LCCC.h"
#include "../IntegerMatrix.h" 
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static void reduce_all_differentials(Komplex *k) {
  if (k == NULL || k->differentials == NULL) return;

  for (int i = 0; i < k->length - 1; i++) {
    if (k->differentials[i] != NULL) {
      CobMatrix_reduce(k->differentials[i]);
    }
  }
}
SmoothingColumn *SmoothingColumn_clone(SmoothingColumn *col) {
  if (col == NULL) return NULL;
  SmoothingColumn *clone = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  if (clone == NULL) return NULL;

  clone->n = col->n;
  if (col->n <= 0) {
    clone->numbers = NULL;
    clone->smoothings = NULL;
    return clone;
  }

  clone->numbers = (int *)malloc((size_t)col->n * sizeof(int));
  clone->smoothings = (Cap **)malloc((size_t)col->n * sizeof(Cap *));
  if (clone->numbers == NULL || clone->smoothings == NULL) {
    free(clone->numbers);
    free(clone->smoothings);
    free(clone);
    return NULL;
  }

  memcpy(clone->numbers, col->numbers, (size_t)col->n * sizeof(int));
  memcpy(clone->smoothings, col->smoothings, (size_t)col->n * sizeof(Cap *));
  return clone;
}

bool SmoothingColumn_equals(SmoothingColumn *a, SmoothingColumn *b) {
  if (a == b) return true;
  if (a == NULL || b == NULL || a->n != b->n) return false;
  for (int i = 0; i < a->n; i++) {
    if (a->numbers[i] != b->numbers[i]) return false;
    if (!Cap_equals(a->smoothings[i], b->smoothings[i])) return false;
  }
  return true;
}

static int quantum_key(const Komplex *k, int h, int idx) {
  if (k == NULL || h < 0 || h >= k->length)
    return 0;

  SmoothingColumn *col = k->chain_groups[h];
  if (col == NULL || col->numbers == NULL || idx < 0 || idx >= col->n)
    return 0;

  return col->numbers[idx];
}

static void free_smoothing_column_owned(SmoothingColumn *col) {
  if (col == NULL)
    return;
  if (col->smoothings != NULL) {
    for (int i = 0; i < col->n; i++) {
      Cap_free(col->smoothings[i]);
    }
  }
  free(col->numbers);
  free(col->smoothings);
  free(col);
}

static void free_komplex_owned(Komplex *k) {
  if (k == NULL)
    return;
  if (k->differentials != NULL) {
    for (int i = 0; i < k->length - 1; i++) {
      CobMatrix_free(k->differentials[i]);
    }
  }
  if (k->chain_groups != NULL) {
    for (int i = 0; i < k->length; i++) {
      free_smoothing_column_owned(k->chain_groups[i]);
    }
  }
  free(k->differentials);
  free(k->chain_groups);
  free(k);
}

static Cap *clone_cap(const Cap *cap) {
  if (cap == NULL)
    return NULL;

  Cap *copy = Cap_create(cap->n, cap->ncycles);
  if (copy == NULL)
    return NULL;
  if (cap->n > 0) {
    memcpy(copy->pairings, cap->pairings, (size_t)cap->n * sizeof(int));
  }
  return copy;
}

static bool komplex_boundary_size(const Komplex *k, int *boundary_size,
                                  char *reason, size_t reason_size) {
  int boundary = -1;
  for (int h = 0; h < k->length; h++) {
    SmoothingColumn *col = k->chain_groups[h];
    if (col == NULL)
      continue;
    for (int i = 0; i < col->n; i++) {
      Cap *cap = col->smoothings[i];
      if (cap == NULL) {
        snprintf(reason, reason_size,
                 "Generator C_%d[%d] has no smoothing object.", h, i);
        return false;
      }
      if (boundary == -1) {
        boundary = cap->n;
      } else if (cap->n != boundary) {
        snprintf(reason, reason_size,
                 "Generator C_%d[%d] has boundary=%d, inconsistent with boundary=%d elsewhere in the complex.",
                 h, i, cap->n, boundary);
        return false;
      }
    }
  }

  *boundary_size = boundary < 0 ? 0 : boundary;
  return true;
}
static Cap *build_braid_closure_cap(int boundary_size) {
  if (boundary_size < 0 || (boundary_size % 2) != 0) return NULL;
  Cap *closure = Cap_create(boundary_size, 0);
  if (closure == NULL) return NULL;

  int half = boundary_size / 2;
  for (int i = 0; i < half; i++) {
    int partner = boundary_size - 1 - i;
    closure->pairings[i] = partner;
    closure->pairings[partner] = i;
  }
  return closure;
}
static Cap *close_cap(const Cap *cap, const Cap *closure_cap, char *reason,
                      size_t reason_size) {
  if (cap == NULL) {
    snprintf(reason, reason_size,
             "Tried to close a missing smoothing object.");
    return NULL;
  }
  if (closure_cap == NULL) {
    snprintf(reason, reason_size,
             "No closure cap is available for braid closure.");
    return NULL;
  }
  if (cap->n != closure_cap->n) {
    snprintf(reason, reason_size,
             "Cannot close a cap with boundary=%d using a closure cap with boundary=%d.",
             cap->n, closure_cap->n);
    return NULL;
  }
  if (cap->n == 0)
    return clone_cap(cap);
  return Cap_compose(cap, 0, closure_cap, 0, cap->n, NULL);
}

static SmoothingColumn *close_smoothing_column(const SmoothingColumn *col,
                                               const Cap *closure_cap,
                                               char *reason,
                                               size_t reason_size) {
  if (col == NULL)
    return NULL;

  SmoothingColumn *closed = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  if (closed == NULL)
    return NULL;

  closed->n = col->n;
  if (col->n <= 0) {
    closed->numbers = NULL;
    closed->smoothings = NULL;
    return closed;
  }

  closed->numbers = (int *)malloc((size_t)col->n * sizeof(int));
  closed->smoothings = (Cap **)malloc((size_t)col->n * sizeof(Cap *));
  if (closed->numbers == NULL || closed->smoothings == NULL) {
    free(closed->numbers);
    free(closed->smoothings);
    free(closed);
    return NULL;
  }

  memcpy(closed->numbers, col->numbers, (size_t)col->n * sizeof(int));
  for (int i = 0; i < col->n; i++) {
    closed->smoothings[i] =
        close_cap(col->smoothings[i], closure_cap, reason, reason_size);
    if (closed->smoothings[i] == NULL) {
      for (int j = 0; j < i; j++) {
        Cap_free(closed->smoothings[j]);
      }
      free(closed->numbers);
      free(closed->smoothings);
      free(closed);
      return NULL;
    }
  }
  return closed;
}

static LCCC *close_lccc(const LCCC *lc, const Cap *closure_cap, int boundary_size,
                        char *reason, size_t reason_size) {
  if (lc == NULL)
    return LCCC_createZero();
  if (boundary_size == 0)
    return LCCC_clone(lc);

  LCCC *result = LCCC_createZero();
  for (LCCCTerm *term = lc->head; term != NULL; term = term->next) {
    if (term->cobordism == NULL) {
      snprintf(reason, reason_size,
               "Encountered a missing cobordism term during braid closure.");
      LCCC_free(result);
      return NULL;
    }

    if (term->cobordism->source == NULL || term->cobordism->target == NULL) {
      snprintf(reason, reason_size,
               "Encountered a malformed cobordism term during braid closure.");
      LCCC_free(result);
      return NULL;
    }
    if (term->cobordism->source->n != boundary_size ||
        term->cobordism->target->n != boundary_size) {
      snprintf(reason, reason_size,
               "Cobordism boundary mismatch during braid closure (source=%d, target=%d, expected=%d).",
               term->cobordism->source->n, term->cobordism->target->n,
               boundary_size);
      LCCC_free(result);
      return NULL;
    }

    CannedCobordism *closure_iso =
        CannedCobordismImpl_isomorphism((Cap *)closure_cap);
    if (closure_iso == NULL) {
      snprintf(reason, reason_size,
               "Failed to build the closure isomorphism cobordism.");
      LCCC_free(result);
      return NULL;
    }

    CannedCobordism *closed_cc = term->cobordism->compose_partial(
        term->cobordism, 0, closure_iso, 0, boundary_size);
    CannedCobordism_free(closure_iso);
    if (closed_cc == NULL) {
      snprintf(reason, reason_size,
               "Horizontal composition failed while closing a cobordism term.");
      LCCC_free(result);
      return NULL;
    }

    LCCC *single = LCCC_createSingle(closed_cc, term->coeff);
    LCCC *sum = LCCC_add(result, single);
    LCCC_free(result);
    LCCC_free(single);
    result = sum;
  }

  return result;
}

static CobMatrix *close_cobmatrix(const CobMatrix *m, SmoothingColumn *source,
                                  SmoothingColumn *target,
                                  const Cap *closure_cap, int boundary_size,
                                  char *reason, size_t reason_size) {
  if (m == NULL)
    return NULL;

  CobMatrix *closed = CobMatrix_create(source, target, true);
  if (closed == NULL)
    return NULL;

  for (int row = 0; row < m->target->n; row++) {
    for (MatrixEntry *entry = m->entries[row].head; entry != NULL;
         entry = entry->next) {
      LCCC *closed_value = close_lccc(entry->value, closure_cap, boundary_size,
                                      reason, reason_size);
      if (closed_value == NULL) {
        CobMatrix_free(closed);
        return NULL;
      }
      CobMatrix_putEntry(closed, row, entry->column_index, closed_value);
    }
  }

  return closed;
}

static Komplex *close_braid_komplex(const Komplex *open, char *reason,
                                    size_t reason_size) {
  int boundary_size = 0;
  if (!komplex_boundary_size(open, &boundary_size, reason, reason_size)) {
    return NULL;
  }

  Cap *closure_cap = build_braid_closure_cap(boundary_size);
  if (closure_cap == NULL) {
    snprintf(reason, reason_size,
             "Failed to build a braid closure cap for boundary size %d.",
             boundary_size);
    return NULL;
  }

  Komplex *closed = Komplex_create(open->length);
  if (closed == NULL) {
    Cap_free(closure_cap);
    return NULL;
  }

  for (int h = 0; h < open->length; h++) {
    closed->chain_groups[h] = close_smoothing_column(
        open->chain_groups[h], closure_cap, reason, reason_size);
    if (open->chain_groups[h] != NULL && closed->chain_groups[h] == NULL) {
      Cap_free(closure_cap);
      free_komplex_owned(closed);
      return NULL;
    }
  }

  for (int i = 0; i < open->length - 1; i++) {
    closed->differentials[i] = close_cobmatrix(
        open->differentials[i], closed->chain_groups[i], closed->chain_groups[i + 1],
        closure_cap, boundary_size, reason, reason_size);
    if (open->differentials[i] != NULL && closed->differentials[i] == NULL) {
      Cap_free(closure_cap);
      free_komplex_owned(closed);
      return NULL;
    }
  }

  Cap_free(closure_cap);
  return closed;
}

static bool komplex_has_only_closed_objects(const Komplex *k, char *reason,
                                            size_t reason_size) {
  for (int h = 0; h < k->length; h++) {
    SmoothingColumn *col = k->chain_groups[h];
    if (col == NULL)
      continue;
    for (int i = 0; i < col->n; i++) {
      Cap *cap = col->smoothings[i];
      if (cap == NULL) {
        snprintf(reason, reason_size,
                 "Generator C_%d[%d] has no smoothing object.", h, i);
        return false;
      }
      if (cap->n != 0 || cap->ncycles != 0) {
        snprintf(reason, reason_size,
                 "Generator C_%d[%d] is still an open tangle object (boundary=%d, cycles=%d).",
                 h, i, cap->n, cap->ncycles);
        return false;
      }
    }
  }
  return true;
}

static bool evaluate_lccc_scalar(LCCC *lc, int64_t *out_value, char *reason,
                                 size_t reason_size) {
  *out_value = 0;
  if (lc == NULL || LCCC_isZero(lc))
    return true;

  LCCC *reduced = LCCC_reduce(lc);
  if (reduced == NULL || LCCC_isZero(reduced)) {
    LCCC_free(reduced);
    return true;
  }

  int64_t total = 0;
  for (LCCCTerm *term = reduced->head; term != NULL; term = term->next) {
    CannedCobordismImplData *impl =
        (CannedCobordismImplData *)term->cobordism->impl_data;

    if (impl == NULL) {
      snprintf(reason, reason_size,
               "Encountered a cobordism term with missing implementation data.");
      LCCC_free(reduced);
      return false;
    }
    if (impl->hpower != 0) {
      snprintf(reason, reason_size,
               "Cannot scalarize a term with nonzero h-power (%d).",
               impl->hpower);
      LCCC_free(reduced);
      return false;
    }
    if (impl->top->n != 0 || impl->bottom->n != 0 || impl->top->ncycles != 0 ||
        impl->bottom->ncycles != 0) {
      snprintf(reason, reason_size,
               "Cannot scalarize a term with non-closed boundary data "
               "(top n=%d cycles=%d, bottom n=%d cycles=%d).",
               impl->top->n, impl->top->ncycles, impl->bottom->n,
               impl->bottom->ncycles);
      LCCC_free(reduced);
      return false;
    }

    //Evaluate each closed connected component using the current Bar-Natan/Khovanov closed-surface rules after reduction.
    int64_t surface_value = 1;

    for (int i = 0; i < impl->ncc; i++) {
      int g = impl->genus[i];
      int d = impl->dots[i];

      if (g == 0 && d == 0) {
        surface_value *= 0;
      } else if (g == 0 && d == 1) {
        surface_value *= 1;
      } else if (g == 1 && d == 0) {
        surface_value *= 2;
      } else {
        surface_value *= 0;
      }
    }

    total += term->coeff * surface_value;
  }

  LCCC_free(reduced);
  *out_value = total;
  return true;
}
static int smith_rank(Mat *m) {
  int diag = m->rows < m->cols ? m->rows : m->cols;
  int rank = 0;
  for (int i = 0; i < diag; i++) {
    if (m->matrix[i][i] != 0)
      rank++;
  }
  return rank;
}
static int run_torsion_regression_test(void) {
  int ok = 1;
  Mat *A = NULL;
  Mat *B = NULL;
  Mat *Acopy = NULL;
  Mat *Bcopy = NULL;

  A = createMat(2, 1);
  B = createMat(1, 2);
  if (A == NULL || B == NULL) { ok = 0; goto cleanup; }

  A->matrix[0][0] = 2;
  A->matrix[1][0] = -2;
  B->matrix[0][0] = 1;
  B->matrix[0][1] = 1;

  if (B->matrix[0][0] * A->matrix[0][0] +
      B->matrix[0][1] * A->matrix[1][0] != 0) {
    ok = 0;
    goto cleanup;
  }

  Bcopy = createMat(1, 2);
  Acopy = createMat(2, 1);
  if (Bcopy == NULL || Acopy == NULL) { ok = 0; goto cleanup; }
  for (int i = 0; i < 1; i++)
    for (int j = 0; j < 2; j++) Bcopy->matrix[i][j] = B->matrix[i][j];
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 1; j++) Acopy->matrix[i][j] = A->matrix[i][j];

  toSmithForm(Bcopy);
  int rank_out = smith_rank(Bcopy);
  toSmithForm(Acopy);
  int rank_in = smith_rank(Acopy);

  int torsion_count = 0;
  int64_t torsion_value = 0;
  int diag = Acopy->rows < Acopy->cols ? Acopy->rows : Acopy->cols;
  for (int i = 0; i < diag; i++) {
    int64_t d = Acopy->matrix[i][i];
    if (d > 1 || d < -1) {
      torsion_value = llabs((long long)d);
      torsion_count++;
    }
  }

  int hom_rank = 2 - rank_out - rank_in;

  if (rank_out != 1 || rank_in != 1 || torsion_count != 1 ||
      torsion_value != 2 || hom_rank != 0)
    ok = 0;

cleanup:
  freeMat(A);
  freeMat(B);
  freeMat(Acopy);
  freeMat(Bcopy);
  return ok;
}

static int q_chain_rank(const Komplex *k, int h, int q) {
  if (k == NULL || h < 0 || h >= k->length ||
      k->chain_groups[h] == NULL)
    return 0;

  int rank = 0;
  for (int i = 0; i < k->chain_groups[h]->n; i++) {
    if (quantum_key(k, h, i) == q)
      rank++;
  }
  return rank;
}

static Mat *build_q_slice(const Komplex *k, Mat *full_diff, int h, int q,
                          char *reason, size_t reason_size) {
  int source_rank = q_chain_rank(k, h, q);
  int target_rank = q_chain_rank(k, h + 1, q);

  if (source_rank <= 0 || target_rank <= 0) {
    snprintf(reason, reason_size,
             "Internal error: requested an empty q-slice for d_%d at q=%d.",
             h, q);
    return NULL;
  }

  Mat *sub = createMat(target_rank, source_rank);
  if (sub == NULL) {
    snprintf(reason, reason_size,
             "Out of memory while building the q-slice of d_%d at q=%d.",
             h, q);
    return NULL;
  }

  int sub_col = 0;
  for (int c = 0; c < k->chain_groups[h]->n; c++) {
    if (quantum_key(k, h, c) != q)
      continue;

    int sub_row = 0;
    for (int r = 0; r < k->chain_groups[h + 1]->n; r++) {
      if (quantum_key(k, h + 1, r) != q)
        continue;

      sub->matrix[sub_row][sub_col] = full_diff->matrix[r][c];
      sub_row++;
    }
    sub_col++;
  }

  return sub;
}

static Mat *convert_CobMatrix_to_IntegerMatrix(int h, CobMatrix *cob_diff,
                                                char *reason,
                                                size_t reason_size,
                                                bool *q_preserving) {
  if (cob_diff == NULL || cob_diff->source == NULL || cob_diff->target == NULL) {
    snprintf(reason, reason_size, "Cannot convert a NULL cobordism matrix.");
    return NULL;
  }

  Mat *m = createMat(cob_diff->target->n, cob_diff->source->n);
  if (m == NULL) {
    snprintf(reason, reason_size, "Failed to allocate the integer matrix.");
    return NULL;
  }

  for (int row = 0; row < cob_diff->target->n; row++) {
    for (MatrixEntry *entry = cob_diff->entries[row].head; entry != NULL;
         entry = entry->next) {
      int col = entry->column_index;
      int64_t value = 0;
      if (!evaluate_lccc_scalar(entry->value, &value, reason, reason_size)) {
        freeMat(m);
        return NULL;
      }

      if (value == 0)
        continue;

      int src_q = 0;
      int tgt_q = 0;

      if (cob_diff->source->numbers != NULL && col >= 0 && col < cob_diff->source->n)
        src_q = cob_diff->source->numbers[col];

      if (cob_diff->target->numbers != NULL && row >= 0 && row < cob_diff->target->n)
        tgt_q = cob_diff->target->numbers[row];

      
      if (src_q != tgt_q) {
        if (q_preserving != NULL)
          *q_preserving = false;

        fprintf(stderr,
                "WARNING: differential d_%d has nonzero scalar entry changing q: "
                "row=%d col=%d value=%lld src_q=%d tgt_q=%d\n",
                h, row, col, (long long)value, src_q, tgt_q);
      }

      m->matrix[row][col] += value;
    }
  }

  return m;
}

static void print_braid_word(const int *crossings, const bool *signs, int length) {
  printf("Braid word:");
  for (int i = 0; i < length; i++) {
    int generator = crossings[i] + 1;
    if (!signs[i]) generator = -generator;
    printf(" %d", generator);
  }
  printf("\n");
}

static void print_scan_poincare(const int *free_ranks, int h_count, int min_q,
                                int max_q, int n_plus, int n_minus,
                                const int *torsion_counts,
                                const int64_t *torsion_values,
                                int max_torsion_per_cell) {
  bool first = true;
  int q_count = max_q - min_q + 1;

  printf("P(q,t) = ");
  for (int h = 0; h < h_count; h++) {
    for (int q = min_q; q <= max_q; q++) {
      int idx = h * q_count + (q - min_q);
      int free_rank = free_ranks[idx];
      int true_h = h - n_minus;
      int true_q = q + n_plus - n_minus;

      if (free_rank > 0) {
        if (!first) printf(" + ");
        first = false;
        if (free_rank != 1) printf("%d", free_rank);
        printf("q^%d", true_q);
        if (true_h != 0) printf("t^%d", true_h);
      }

      for (int i = 0; i < torsion_counts[idx]; i++) {
        if (!first) printf(" + ");
        first = false;
        printf("Z_%lld",
               (long long)torsion_values[idx * max_torsion_per_cell + i]);
        printf("q^%d", true_q);
        if (true_h != 0) printf("t^%d", true_h);
      }
    }
  }

  if (first) printf("0");
  printf("\n");
}

static Komplex *build_scan_komplex(int n_strands, const int *crossings,
                                   const bool *signs, int length) {
  Komplex *current = Komplex_identityBraid(n_strands);
  if (current == NULL) return NULL;

  for (int i = 0; i < length; i++) {
    Komplex *cross = Komplex_singleCrossing(n_strands, crossings[i], signs[i]);
    Komplex *next = Komplex_compose_tangles(current, cross, n_strands);
    Komplex_free(current);
    Komplex_free(cross);
    if (next == NULL) return NULL;

    current = next;

    Komplex_deloop(current);

    Komplex_greedyReduce(current);
  }

  return current;
}

int main(int argc, char **argv) {
  if (!run_torsion_regression_test()) {
    fprintf(stderr, "torsion regression test failed\n");
    return 3;
  }
  int n_strands = 0;
  int *crossings = NULL;
  bool *signs = NULL;
  int length = 0;
  int strand_count = 2;
  int start = 1;
  while (start < argc) {
    if (strcmp(argv[start], "--quiet") == 0) {
      start++;
      continue;
    }

    if (strcmp(argv[start], "--n-strands") == 0) {
      if (start + 1 >= argc) return 2;
      strand_count = atoi(argv[start + 1]);
      start += 2;
      continue;
    }

    break;
  }

  n_strands = strand_count;

  if (start >= argc) {
    length = 3;
    n_strands = 2;
    crossings = (int *)malloc(3 * sizeof(int));
    signs = (bool *)malloc(3 * sizeof(bool));
    if (crossings == NULL || signs == NULL) {
      free(crossings);
      free(signs);
      return 1;
    }

    for (int i = 0; i < 3; i++) {
      crossings[i] = 0;
      signs[i] = false;
    }
  } else {
    length = argc - start;
    crossings = (int *)malloc((size_t)length * sizeof(int));
    signs = (bool *)malloc((size_t)length * sizeof(bool));
    if (crossings == NULL || signs == NULL) {
      free(crossings);
      free(signs);
      return 1;
    }

    for (int i = 0; i < length; i++) {
      int generator = atoi(argv[start + i]);
      int abs_generator = generator > 0 ? generator : -generator;

      if (generator == 0 || abs_generator >= n_strands) {
        fprintf(stderr,
                "Invalid braid generator %d for %d strands.\n",
                generator, n_strands);
        free(crossings);
        free(signs);
        return 2;
      }

      crossings[i] = abs_generator - 1;
      signs[i] = generator > 0;
    }
  }

  int n_plus = 0;
  int n_minus = 0;
  for (int i = 0; i < length; i++) {
    if (signs[i]) n_plus++;
    else n_minus++;
  }

  printf("Scanning Khovanov harness\n");
  print_braid_word(crossings, signs, length);

  Komplex *result = build_scan_komplex(n_strands, crossings, signs, length);
  free(crossings);
  free(signs);

  if (result == NULL) return 1;

  char reason[256];

  Komplex *closed = close_braid_komplex(result, reason, sizeof(reason));
  if (closed == NULL) {
    printf("Blocked: %s\n", reason);
    free_komplex_owned(result);
    return 1;
  }

  Komplex_deloop(closed);
  reduce_all_differentials(closed);

  Komplex_greedyReduce(closed);
  reduce_all_differentials(closed);

  if (!komplex_has_only_closed_objects(closed, reason, sizeof(reason))) {
    printf("Blocked: %s\n", reason);
    free_komplex_owned(closed);
    free_komplex_owned(result);
    return 1;
  }

  Mat **full_diffs = NULL;
  if (closed->length > 1) {
    full_diffs = (Mat **)calloc((size_t)(closed->length - 1), sizeof(Mat *));
    if (full_diffs == NULL) {
      printf("Blocked: out of memory while allocating integer differentials.\n");
      free_komplex_owned(closed);
      free_komplex_owned(result);
      return 1;
    }
  }

  bool ok = true;
  bool q_preserving = true;
  for (int i = 0; i < closed->length - 1; i++) {
    full_diffs[i] = convert_CobMatrix_to_IntegerMatrix(
        i, closed->differentials[i], reason, sizeof(reason), &q_preserving);
    if (full_diffs[i] == NULL) {
      ok = false;
      break;
    }
  }

  if (!ok) {
    printf("Blocked: %s\n", reason);
  } else if (!q_preserving) {
    printf("Blocked: closed scalar complex is not q-preserving under the current grading labels.\n");
    printf("Skipping bigraded homology/Poincare extraction because same-q submatrices are not subcomplexes.\n");
  } else {
    int min_q = 10000;
    int max_q = -10000;
    int max_torsion_per_cell = 1;

    for (int h = 0; h < closed->length; h++) {
      if (closed->chain_groups[h] != NULL &&
          closed->chain_groups[h]->n > max_torsion_per_cell)
        max_torsion_per_cell = closed->chain_groups[h]->n;
      if (!closed->chain_groups[h]) continue;
      for (int i = 0; i < closed->chain_groups[h]->n; i++) {
        int q = quantum_key(closed, h, i);
        if (q < min_q) min_q = q;
        if (q > max_q) max_q = q;
      }
    }

    if (max_q < min_q) {
      printf("P(q,t) = 0\n");
    } else {
      int q_count = max_q - min_q + 1;
      int cell_count = closed->length * q_count;
      int *free_ranks = (int *)calloc((size_t)cell_count, sizeof(int));
      int *torsion_counts = (int *)calloc((size_t)cell_count, sizeof(int));
      int64_t *poincare_torsion_values =
          (int64_t *)calloc(
              (size_t)cell_count * (size_t)max_torsion_per_cell,
              sizeof(int64_t));

      if (free_ranks == NULL || torsion_counts == NULL ||
          poincare_torsion_values == NULL) {
        printf("Blocked: out of memory while building the Poincare summary.\n");
        free(free_ranks);
        free(torsion_counts);
        free(poincare_torsion_values);

        if (full_diffs != NULL) {
          for (int i = 0; i < closed->length - 1; i++) freeMat(full_diffs[i]);
          free(full_diffs);
        }

        free_komplex_owned(closed);
        free_komplex_owned(result);
        return 1;
      }

      for (int q = min_q; q <= max_q; q++) {
        bool has_q = false;
        for (int h = 0; h < closed->length; h++) {
          if (!closed->chain_groups[h]) continue;
          for (int i = 0; i < closed->chain_groups[h]->n; i++) {
            if (quantum_key(closed, h, i) == q) has_q = true;
          }
        }

        if (!has_q) continue;

        for (int h = 0; h < closed->length; h++) {
          int chain_rank = q_chain_rank(closed, h, q);
          if (chain_rank == 0)
            continue;

          int source_rank = q_chain_rank(closed, h - 1, q);
          int target_rank = q_chain_rank(closed, h + 1, q);
          int cell_idx = h * q_count + (q - min_q);

          Mat *incoming = NULL;
          Mat *outgoing = NULL;

          if (source_rank > 0) {
            incoming = build_q_slice(closed, full_diffs[h - 1], h - 1, q,
                                     reason, sizeof(reason));
            if (incoming == NULL) {
              printf("Blocked: %s\n", reason);
              ok = false;
              break;
            }
          }

          if (target_rank > 0) {
            outgoing = build_q_slice(closed, full_diffs[h], h, q,
                                     reason, sizeof(reason));
            if (outgoing == NULL) {
              printf("Blocked: %s\n", reason);
              freeMat(incoming);
              ok = false;
              break;
            }
          }

          int rank_out = 0;
          int rank_in = 0;
          int torsion_count = 0;

          if (outgoing != NULL) {
            toSmithForm(outgoing);
            if (!isDiag(outgoing)) {
              printf("Blocked: failed to diagonalize outgoing q-slice at h=%d, q=%d.\n",
                     h, q);
              freeMat(outgoing);
              freeMat(incoming);
              ok = false;
              break;
            }
            rank_out = smith_rank(outgoing);
          }

          if (incoming != NULL) {
            toSmithForm(incoming);
            if (!isDiag(incoming)) {
              printf("Blocked: failed to diagonalize incoming q-slice at h=%d, q=%d.\n",
                     h, q);
              freeMat(outgoing);
              freeMat(incoming);
              ok = false;
              break;
            }

            rank_in = smith_rank(incoming);
            int diag = incoming->rows < incoming->cols
                           ? incoming->rows
                           : incoming->cols;

            for (int i = 0; i < diag; i++) {
              int64_t value = incoming->matrix[i][i];
              if (value < 0) value = -value;

              if (value > 1) {
                if (torsion_count >= max_torsion_per_cell) {
                  printf("Blocked: torsion-storage bound exceeded at h=%d, q=%d.\n",
                         h, q);
                  ok = false;
                  break;
                }

                poincare_torsion_values[
                    cell_idx * max_torsion_per_cell + torsion_count] = value;
                torsion_count++;
              }
            }

            if (!ok) {
              freeMat(outgoing);
              freeMat(incoming);
              break;
            }
          }

          int hom_rank = chain_rank - rank_out - rank_in;
          if (hom_rank < 0) {
            printf("Blocked: computed a negative homology rank at h=%d, q=%d.\n",
                   h, q);
            freeMat(outgoing);
            freeMat(incoming);
            ok = false;
            break;
          }

          /* THE TRUE KHOVANOV SHIFT */
          int true_h = h - n_minus;
          int true_q = q + n_plus - n_minus;

          free_ranks[cell_idx] = hom_rank;
          torsion_counts[cell_idx] = torsion_count;

          if (hom_rank > 0) {
            printf("rank H^{%d, %d} = %d\n", true_h, true_q, hom_rank);
          } else if (torsion_count == 0) {
            printf("rank H^{%d, %d} = 0\n", true_h, true_q);
          }

          for (int i = 0; i < torsion_count; i++) {
            int64_t value = poincare_torsion_values[
                cell_idx * max_torsion_per_cell + i];
            printf("torsion H^{%d, %d} = Z_%lld\n",
                   true_h, true_q, (long long)value);
          }

          freeMat(outgoing);
          freeMat(incoming);
        }

        if (!ok) break;
      }

      if (ok) {
        printf("\n--- Poincare Polynomial ---\n");
        print_scan_poincare(free_ranks, closed->length, min_q, max_q, n_plus,
                            n_minus, torsion_counts, poincare_torsion_values,
                            max_torsion_per_cell);
      }

      free(free_ranks);
      free(torsion_counts);
      free(poincare_torsion_values);
    }
  }

  if (full_diffs != NULL) {
    for (int i = 0; i < closed->length - 1; i++) freeMat(full_diffs[i]);
    free(full_diffs);
  }

  free_komplex_owned(closed);
  free_komplex_owned(result);
  return ok ? 0 : 1;
}