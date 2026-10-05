/* Standalone FastKh scanner. Default stdout is one free-rank/field polynomial.
 * Allocation and arithmetic failures terminate before a result is printed.
 * Cap/cobordism ownership remains the supplied pipeline's shared-pointer model. */
#include "../KhMemory.h"
#include "TangleKomplex.h"
#include "PDScanner.h"
#include "DTtoPD.h"
#include "Komplex.h"
#include "Cap.h"
#include "LCCC.h"
#include "../IntegerMatrix.h"
#include "../polynomial/BivariatePoly.h"
#include <stdbool.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
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
static Cap *clone_cap(const Cap *cap) {
  if (cap == NULL)
    return NULL;
  Cap *copy = Cap_create(cap->n, cap->ncycles);

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
  SmoothingColumn *closed = SmoothingColumn_create();

  closed->n = col->n;
  if (col->n <= 0) {
    closed->numbers = NULL;
    closed->smoothings = NULL;
    return closed;
  }
  closed->numbers = (int *)kh_malloc((size_t)col->n * sizeof(int));
  closed->smoothings = (Cap **)kh_malloc((size_t)col->n * sizeof(Cap *));

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
    CannedCobordism_free(closed_cc);
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

  for (int h = 0; h < open->length; h++) {
    closed->chain_groups[h] = close_smoothing_column(
      open->chain_groups[h], closure_cap, reason, reason_size);
    if (open->chain_groups[h] != NULL && closed->chain_groups[h] == NULL) {
      Cap_free(closure_cap);
      Komplex_free(closed);
      return NULL;
    }
  }
  for (int i = 0; i < open->length - 1; i++) {
    closed->differentials[i] = close_cobmatrix(
      open->differentials[i], closed->chain_groups[i], closed->chain_groups[i + 1],
      closure_cap, boundary_size, reason, reason_size);
    if (open->differentials[i] != NULL && closed->differentials[i] == NULL) {
      Cap_free(closure_cap);
      Komplex_free(closed);
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
static int64_t normalize_mod_i64(int64_t value, int p) {
  int64_t r = value % (int64_t)p;
  if (r < 0) r += p;
  return r;
}
static bool evaluate_lccc_scalar(LCCC *lc, int64_t *out_value, char *reason,
  size_t reason_size) {
  *out_value = 0;
  if (lc == NULL || LCCC_isZero(lc))
    return true;
  LCCC *reduced = LCCC_reduce(lc);

  if (LCCC_isZero(reduced)) {
    LCCC_free(reduced);
    return true;
  }
  bool field_mode = LCCC_getCoefficientMode() == KH_COEFF_FP;
  int prime = LCCC_getCoefficientModulus();
  int64_t total_z = 0;
  int total_fp = 0;
  for (LCCCTerm *term = reduced->head; term != NULL; term = term->next) {
    if (term->cobordism == NULL) kh_fatal("missing scalar cobordism");
    CannedCobordismImplData *impl =
      (CannedCobordismImplData *)term->cobordism->impl_data;
    if (impl == NULL || impl->top == NULL || impl->bottom == NULL) {
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
    if (field_mode) {
      int term_value = LCCC_coeffNormalize(term->coeff);
      for (int i = 0; i < impl->ncc; i++) {
        int g = impl->genus[i];
        int d = impl->dots[i];
        if (g < 0 || d < 0) kh_fatal("negative genus or dot count");
        int surface_factor = 0;
        if (g == 0 && d == 1)
          surface_factor = 1;
        else if (g == 1 && d == 0)
          surface_factor = 2;
        term_value = LCCC_coeffMultiply(term_value, surface_factor);
        if (LCCC_coeffIsZero(term_value)) break;
      }
      total_fp = LCCC_coeffAdd(total_fp, term_value);
    } else {
      /* Preserve the existing integral scalarization path. */
      int64_t surface_value = 1;
      for (int i = 0; i < impl->ncc; i++) {
        int g = impl->genus[i];
        int d = impl->dots[i];
        if (g < 0 || d < 0) kh_fatal("negative genus or dot count");
        if (g == 0 && d == 0) {
          surface_value *= 0;
        } else if (g == 0 && d == 1) {
          surface_value *= 1;
        } else if (g == 1 && d == 0) {
          surface_value = kh_mul64(surface_value, 2);
        } else {
          surface_value *= 0;
        }
      }
      total_z = kh_add64(total_z, kh_mul64(term->coeff, surface_value));
    }
  }
  LCCC_free(reduced);
  *out_value = field_mode ? normalize_mod_i64(total_fp, prime) : total_z;
  return true;
}

/* Finite-field elimination consumes the q-block; no second matrix copy. */
static int rank_mod_prime(Mat *m, int prime) {
  int rank = 0;
  for (int col = 0; col < m->cols && rank < m->rows; ++col) {
    int pivot = rank;
    while (pivot < m->rows && m->matrix[pivot][col] == 0) ++pivot;
    if (pivot == m->rows) continue;
    int64_t *tmp = m->matrix[rank];
    m->matrix[rank] = m->matrix[pivot];
    m->matrix[pivot] = tmp;
    int inverse;
    if (!LCCC_coeffInverse((int)m->matrix[rank][col], &inverse))
      kh_fatal("noninvertible finite-field pivot");
    for (int c = col; c < m->cols; ++c)
      m->matrix[rank][c] = normalize_mod_i64(m->matrix[rank][c] * inverse, prime);
    for (int r = rank + 1; r < m->rows; ++r) {
      int64_t factor = m->matrix[r][col];
      if (!factor) continue;
      for (int c = col; c < m->cols; ++c)
        m->matrix[r][c] = normalize_mod_i64(
            m->matrix[r][c] - factor * m->matrix[rank][c], prime);
    }
    ++rank;
  }
  return rank;
}

typedef struct {
  int count;
  int64_t *values;
} TorsionList;

typedef struct {
  int length, nq;
  int *qs;
  int **q_index, **local_index;
  int *sizes, *ranks;
  Mat **blocks;
  TorsionList *torsion; /* indexed by differential; belongs to next degree */
} ScalarComplex;

static int compare_ints(const void *a, const void *b) {
  int x = *(const int *)a, y = *(const int *)b;
  return (x > y) - (x < y);
}
static int find_q(const int *qs, int n, int q) {
  int lo = 0, hi = n;
  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    if (qs[mid] < q) lo = mid + 1;
    else hi = mid;
  }
  if (lo == n || qs[lo] != q) kh_fatal("missing quantum block");
  return lo;
}
static size_t cell(const ScalarComplex *s, int h, int qi) {
  return (size_t)h * (size_t)s->nq + (size_t)qi;
}
static void scalar_free(ScalarComplex *s) {
  for (int h = 0; h < s->length; ++h) {
    free(s->q_index[h]);
    free(s->local_index[h]);
    for (int qi = 0; qi < s->nq; ++qi) {
      size_t i = cell(s, h, qi);
      freeMat(s->blocks[i]);
      if (s->torsion) free(s->torsion[i].values);
    }
  }
  free(s->qs); free(s->q_index); free(s->local_index);
  free(s->sizes); free(s->ranks); free(s->blocks); free(s->torsion);
}

/* Each sparse cobordism entry is evaluated once and routed directly to its
 * q-block. Null blocks represent zero maps with dimensions stored in sizes. */
static ScalarComplex scalar_build(Komplex *k, bool keep_torsion) {
  if (!k || k->length < 1) kh_fatal("invalid complex length");
  ScalarComplex s = {0};
  s.length = k->length;
  size_t total = 0;
  for (int h = 0; h < k->length; ++h) {
    SmoothingColumn *c = k->chain_groups[h];
    if (!c) continue;
    if (c->n < 0 || (c->n && (!c->numbers || !c->smoothings)))
      kh_fatal("invalid smoothing column");
    if ((size_t)c->n > (size_t)INT_MAX - total)
      kh_fatal("too many generators for quantum indexing");
    total += (size_t)c->n;
  }
  s.qs = kh_malloc(kh_size_mul(total, sizeof(int)));
  size_t at = 0;
  for (int h = 0; h < k->length; ++h) {
    SmoothingColumn *c = k->chain_groups[h];
    if (c) for (int j = 0; j < c->n; ++j) s.qs[at++] = c->numbers[j];
  }
  qsort(s.qs, total, sizeof(int), compare_ints);
  for (size_t i = 0; i < total; ++i)
    if (!s.nq || s.qs[i] != s.qs[s.nq - 1]) s.qs[s.nq++] = s.qs[i];
  size_t cells = kh_size_mul((size_t)s.length, (size_t)s.nq);
  s.sizes = kh_calloc(cells, sizeof(int));
  s.ranks = kh_calloc(cells, sizeof(int));
  s.blocks = kh_calloc(cells, sizeof(Mat *));
  if (keep_torsion) s.torsion = kh_calloc(cells, sizeof(TorsionList));
  s.q_index = kh_calloc((size_t)s.length, sizeof(int *));
  s.local_index = kh_calloc((size_t)s.length, sizeof(int *));
  for (int h = 0; h < s.length; ++h) {
    SmoothingColumn *c = k->chain_groups[h];
    if (!c) continue;
    s.q_index[h] = kh_calloc((size_t)c->n, sizeof(int));
    s.local_index[h] = kh_calloc((size_t)c->n, sizeof(int));
    for (int j = 0; j < c->n; ++j) {
      int qi = find_q(s.qs, s.nq, c->numbers[j]);
      s.q_index[h][j] = qi;
      s.local_index[h][j] = s.sizes[cell(&s, h, qi)]++;
    }
  }
  const bool fp = LCCC_getCoefficientMode() == KH_COEFF_FP;
  const int prime = LCCC_getCoefficientModulus();
  for (int h = 0; h < s.length - 1; ++h) {
    CobMatrix *d = k->differentials[h];
    SmoothingColumn *src = k->chain_groups[h], *dst = k->chain_groups[h + 1];
    if (!d) {
      if (src && dst && src->n && dst->n)
        kh_fatal("missing differential between nonempty chain groups");
      continue;
    }
    if (!SmoothingColumn_equals(d->source, src) ||
        !SmoothingColumn_equals(d->target, dst))
      kh_fatal("differential endpoints disagree with chain groups");
    if (!src || !dst || (dst->n && !d->entries))
      kh_fatal("malformed differential");
    for (int r = 0; r < dst->n; ++r) {
      for (MatrixEntry *e = d->entries[r].head; e; e = e->next) {
        int c = e->column_index;
        if (c < 0 || c >= src->n || !e->value)
          kh_fatal("invalid sparse differential entry");
        int64_t value;
        char reason[256] = {0};
        if (!evaluate_lccc_scalar(e->value, &value, reason, sizeof(reason)))
          kh_fatal(reason[0] ? reason : "scalar evaluation failed");
        if (!value) continue;
        int qi = s.q_index[h][c];
        if (qi != s.q_index[h + 1][r])
          kh_fatal("nonzero scalar differential changes quantum grading");
        size_t idx = cell(&s, h, qi);
        if (!s.blocks[idx]) {
          s.blocks[idx] = createMat(s.sizes[cell(&s, h + 1, qi)], s.sizes[idx]);
          if (!s.blocks[idx]) kh_fatal("could not create quantum block");
        }
        int64_t *entry = &s.blocks[idx]->matrix[s.local_index[h + 1][r]][s.local_index[h][c]];
        *entry = fp ? normalize_mod_i64(*entry + value, prime) : kh_add64(*entry, value);
      }
    }
  }
  return s;
}

/* Exact check in the selected coefficient system, before destructive reduction.
 * Fixed-width integral overflow is an explicit failure, never a passed check. */
static void scalar_verify(const ScalarComplex *s) {
  bool fp = LCCC_getCoefficientMode() == KH_COEFF_FP;
  int prime = LCCC_getCoefficientModulus();
  for (int qi = 0; qi < s->nq; ++qi) {
    for (int h = 0; h < s->length - 2; ++h) {
      Mat *a = s->blocks[cell(s, h, qi)], *b = s->blocks[cell(s, h + 1, qi)];
      if (!a || !b) continue;
      if (b->cols != a->rows) kh_fatal("incompatible consecutive scalar blocks");
      for (int r = 0; r < b->rows; ++r) {
        for (int c = 0; c < a->cols; ++c) {
          int64_t sum = 0;
          for (int j = 0; j < a->rows; ++j) {
            if (!b->matrix[r][j] || !a->matrix[j][c]) continue;
            int64_t product = kh_mul64(b->matrix[r][j], a->matrix[j][c]);
            sum = fp ? normalize_mod_i64(sum + product, prime) : kh_add64(sum, product);
          }
          if (sum) kh_fatal("scalar differential does not satisfy d^2 = 0");
        }
      }
    }
  }
}

static BivariatePoly *scalar_homology(ScalarComplex *s, int n_plus, int n_minus) {
  bool fp = LCCC_getCoefficientMode() == KH_COEFF_FP;
  int prime = LCCC_getCoefficientModulus();
  BivariatePoly *poly = bp_create();

  for (int qi = 0; qi < s->nq; ++qi) {
    for (int h = 0; h < s->length - 1; ++h) {
      size_t idx = cell(s, h, qi);
      Mat *m = s->blocks[idx];
      if (!m) continue;
      if (fp) s->ranks[idx] = rank_mod_prime(m, prime);
      else {
        toSmithForm(m);
        if (!isDiag(m)) kh_fatal("integer differential reduction is not diagonal");
        int n = m->rows < m->cols ? m->rows : m->cols;
        int nt = 0;
        for (int j = 0; j < n; ++j) {
          int64_t v = m->matrix[j][j];
          if (v == INT64_MIN) kh_fatal("torsion order exceeds int64_t");
          if (v) ++s->ranks[idx];
          if (v > 1 || v < -1) ++nt;
        }
        if (s->torsion && nt) {
          TorsionList *t = &s->torsion[idx];
          t->values = kh_calloc((size_t)nt, sizeof(int64_t));
          for (int j = 0; j < n; ++j) {
            int64_t v = m->matrix[j][j];
            if (v < 0) v = -v;
            if (v > 1) t->values[t->count++] = v;
          }
        }
      }
      freeMat(m);
      s->blocks[idx] = NULL;
    }
    for (int h = 0; h < s->length; ++h) {
      size_t idx = cell(s, h, qi);
      int rank_in = h ? s->ranks[cell(s, h - 1, qi)] : 0;
      int64_t betti = (int64_t)s->sizes[idx] - rank_in - s->ranks[idx];
      if (betti < 0) kh_fatal("negative homology rank");
      bp_add_term(poly, kh_int((int64_t)s->qs[qi] + n_plus - n_minus),
                  kh_int((int64_t)h - n_minus), kh_int(betti));
    }
  }
  return poly;
}

static void print_result(const BivariatePoly *poly, const ScalarComplex *s,
                         int n_plus, int n_minus, bool detailed) {
  bool fp = LCCC_getCoefficientMode() == KH_COEFF_FP;
  if (detailed) {
    if (fp) printf("Kh_F%d(q,t) = ", LCCC_getCoefficientModulus());
    else printf("Kh_free(q,t) = ");
  }
  bp_print(poly);
  if (detailed) {
    for (int i = 0; i < poly->num_terms; ++i) {
      const Monomial *m = &poly->terms[i];
      if (fp) printf("dim_F%d H^{%d, %d} = %d\n", LCCC_getCoefficientModulus(), m->t_exp, m->q_exp, m->coeff);
      else printf("rank H^{%d, %d} = %d\n", m->t_exp, m->q_exp, m->coeff);
    }
    if (!fp) {
      bool any = false;
      for (int h = 0; h < s->length - 1; ++h)
        for (int qi = 0; qi < s->nq; ++qi) {
          const TorsionList *t = &s->torsion[cell(s, h, qi)];
          for (int j = 0; j < t->count; ++j) {
            printf("torsion H^{%d, %d} = Z_%lld\n", kh_int((int64_t)h + 1 - n_minus),
                   kh_int((int64_t)s->qs[qi] + n_plus - n_minus), (long long)t->values[j]);
            any = true;
          }
        }
      if (!any) puts("Torsion terms: none");
    }
  }
  if (fflush(stdout) == EOF || ferror(stdout)) kh_fatal("could not write result");
}

static Komplex *build_scan_komplex(int strands, const int *generators, int count,
                                  bool verify) {
  Komplex *current = Komplex_identityBraid(strands);
  if (!current) kh_fatal("could not construct identity braid");
  for (int i = 0; i < count; ++i) {
    int g = generators[i];
    int crossing = (g < 0 ? -g : g) - 1;
    Komplex *cross = Komplex_singleCrossing(strands, crossing, g > 0);
    if (!cross) kh_fatal("could not construct braid crossing");
    Komplex *next = Komplex_compose_tangles(current, cross, strands);
    if (!next) kh_fatal("could not compose braid crossing");
    Komplex_free(current); Komplex_free(cross);
    current = next;
    Komplex_deloop(current);
    reduce_all_differentials(current);
    Komplex_greedyReduce(current);
    reduce_all_differentials(current);
    if (verify && !Komplex_verify_d_squared(current)) kh_fatal("braid d^2 check failed");
  }
  return current;
}

static void usage(FILE *out) {
  fputs("Usage: FullScanning [options] --dt CODE\n"
        "       FullScanning [options] --pd 'PD[...]' [--signs LIST | --javakh-signs]\n"
        "       FullScanning [options] --n-strands N [signed braid generators]\n"
        "Options:\n"
        "  --coeff Z|F<p>     Integral coefficients (default), or prime field\n"
        "  --mod p           Alias for --coeff F<p>\n"
        "  --format polynomial|detailed  Default: one polynomial on stdout\n"
        "  --quiet           Alias for --format polynomial; errors remain on stderr\n"
        "  --no-reorder      Disable PD/DT crossing reordering\n"
        "  --no-d2-check     Skip chain-condition checks (grading checks remain)\n"
        "  --help            Show this help\n"
        "Braid generators default to two strands if --n-strands is omitted.\n"
        "An explicit strand count with no generators denotes the identity braid.\n", out);
}
static int parse_int(const char *text, const char *what) {
  errno = 0;
  char *end = NULL;
  long value = strtol(text, &end, 10);
  if (errno == ERANGE || end == text || *end || value < INT_MIN || value > INT_MAX) {
    fprintf(stderr, "Invalid %s: %s\n", what, text);
    exit(2);
  }
  return (int)value;
}
static const char *option_value(int argc, char **argv, int *i) {
  if (*i + 1 >= argc) {
    fprintf(stderr, "%s requires a value.\n", argv[*i]);
    exit(2);
  }
  return argv[++*i];
}
static _Noreturn void input_error(const char *message) {
  fprintf(stderr, "%s\n", message);
  exit(2);
}

int main(int argc, char **argv) {
  LCCC_setCoefficientIntegers();
  const char *pd = NULL, *dt = NULL, *signs = NULL;
  bool java_signs = false, reorder = true, verify = true, detailed = false;
  bool strands_set = false, reorder_set = false;
  int strands = 2, count = 0, n_plus = 0, n_minus = 0;
  int *generators = kh_calloc((size_t)argc, sizeof(int));
  for (int i = 1; i < argc; ++i) {
    const char *arg = argv[i];
    if (!strcmp(arg, "--help")) { usage(stdout); free(generators); return 0; }
    else if (!strcmp(arg, "--quiet")) detailed = false;
    else if (!strcmp(arg, "--format")) {
      const char *v = option_value(argc, argv, &i);
      if (!strcmp(v, "detailed")) detailed = true;
      else if (!strcmp(v, "polynomial")) detailed = false;
      else input_error("--format requires polynomial or detailed.");
    } else if (!strcmp(arg, "--coeff") || !strcmp(arg, "--mod")) {
      const char *v = option_value(argc, argv, &i);
      if (!strcmp(arg, "--coeff") && (!strcmp(v, "Z") || !strcmp(v, "z")))
        LCCC_setCoefficientIntegers();
      else {
        if (!strcmp(arg, "--coeff")) {
          if (*v != 'F' && *v != 'f') input_error("--coeff requires Z or F<p>.");
          ++v;
        }
        int prime = parse_int(v, "prime modulus");
        if (!LCCC_setCoefficientModPrime(prime)) input_error("The modulus must be prime.");
      }
    } else if (!strcmp(arg, "--n-strands")) {
      strands = parse_int(option_value(argc, argv, &i), "strand count");
      if (strands < 1 || strands > INT_MAX / 2) input_error("Strand count is out of range.");
      strands_set = true;
    } else if (!strcmp(arg, "--pd")) {
      if (pd) input_error("Repeated --pd input.");
      pd = option_value(argc, argv, &i);
    } else if (!strcmp(arg, "--dt")) {
      if (dt) input_error("Repeated --dt input.");
      dt = option_value(argc, argv, &i);
    } else if (!strcmp(arg, "--signs")) signs = option_value(argc, argv, &i);
    else if (!strcmp(arg, "--javakh-signs")) java_signs = true;
    else if (!strcmp(arg, "--no-reorder")) { reorder = false; reorder_set = true; }
    else if (!strcmp(arg, "--no-d2-check")) verify = false;
    else {
      if (!strncmp(arg, "--", 2)) input_error("Unknown option; use --help.");
      generators[count++] = parse_int(arg, "braid generator");
    }
  }
  if (pd && dt) input_error("--pd and --dt are mutually exclusive.");
  if ((pd || dt) && (count || strands_set)) input_error("PD/DT input cannot be combined with braid input.");
  if (!pd && (signs || java_signs)) input_error("--signs and --javakh-signs require --pd.");
  if (signs && java_signs) input_error("Choose --signs or --javakh-signs, not both.");
  if (!pd && !dt && reorder_set) input_error("--no-reorder requires PD or DT input.");
  if (!pd && !dt && !count && !strands_set) {
    usage(stderr); free(generators); return 2;
  }
  char reason[256] = {0};
  Komplex *open = NULL, *closed = NULL;
  if (pd || dt) {
    PDDiagram diagram = {0};
    bool parsed = dt ? DTtoPD_fromAlphabetical(&diagram, dt, reason, sizeof(reason))
                     : PDDiagram_parse(&diagram, pd, signs, java_signs, reason, sizeof(reason));
    if (!parsed) input_error(reason[0] ? reason : "invalid knot diagram");
    for (int i = 0; i < diagram.crossing_count; ++i) {
      if (diagram.signs[i] == 1) ++n_plus;
      else if (diagram.signs[i] == -1) ++n_minus;
      else input_error("Crossing signs must be +1 or -1.");
    }
    closed = PDScanner_build(&diagram, reorder, verify, reason, sizeof(reason));
    PDDiagram_free(&diagram);
    if (!closed) kh_fatal(reason[0] ? reason : "PD scan failed");
    /* PDScanner already closes and simplifies its last crossing. */
  } else {
    for (int i = 0; i < count; ++i) {
      int g = generators[i];
      if (!g || g == INT_MIN || g >= strands || g <= -strands)
        input_error("Braid generator is out of range for the selected strand count.");
      if (g > 0) ++n_plus; else ++n_minus;
    }
    open = build_scan_komplex(strands, generators, count, verify);
    closed = close_braid_komplex(open, reason, sizeof(reason));
    if (!closed) kh_fatal(reason[0] ? reason : "braid closure failed");
    Komplex_deloop(closed);
    reduce_all_differentials(closed);
    Komplex_greedyReduce(closed);
    reduce_all_differentials(closed);
  }
  free(generators);
  if (!komplex_has_only_closed_objects(closed, reason, sizeof(reason))) kh_fatal(reason);
  ScalarComplex scalar = scalar_build(closed, detailed && LCCC_getCoefficientMode() == KH_COEFF_Z);
  if (verify) scalar_verify(&scalar);
  BivariatePoly *poly = scalar_homology(&scalar, n_plus, n_minus);
  print_result(poly, &scalar, n_plus, n_minus, detailed);
  bp_free(poly);
  scalar_free(&scalar);
  Komplex_free(closed);
  Komplex_free(open);
  return 0;
}
