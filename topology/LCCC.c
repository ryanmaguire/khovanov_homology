/*
 *  LCCC.c
 *
 *  Linear Combination of Canned Cobordisms with coefficients in Z or F_p.
 *
 *  Purpose:
 *      Implements the coefficient type used in the cobordism-valued
 *      Khovanov complex. Coefficients remain stored as int values, while
 *      arithmetic is interpreted either over Z or over a selected prime
 *      field F_p. Reduction applies the Bar-Natan surface relations in
 *      that same coefficient system.
 *
 *  Compile (standalone test):
 *      gcc -Wall -Wextra -g -c LCCC.c -I.
 */

#include "LCCC.h"
#include "CannedCobordismImpl.h"
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 *  Coefficient mode
 * ================================================================ */

/*
 * One FullScanning run uses one coefficient system. Keeping the selection
 * here avoids changing every existing LCCC/CobMatrix/Komplex function
 * signature. Integer coefficients remain the default, preserving the old
 * behavior unless field mode is explicitly selected.
 */
static KHCoefficientMode g_coeff_mode = KH_COEFF_Z;
static int g_coeff_modulus = 0;

static bool is_prime_int(int p) {
  if (p < 2)
    return false;
  if (p == 2)
    return true;
  if ((p % 2) == 0)
    return false;

  for (int d = 3; d <= p / d; d += 2) {
    if ((p % d) == 0)
      return false;
  }
  return true;
}

static int normalize_i64(int64_t value) {
  if (g_coeff_mode == KH_COEFF_Z) {
    if (value > INT_MAX || value < INT_MIN)
      abort();
    return (int)value;
  }

  int64_t p = (int64_t)g_coeff_modulus;
  int64_t r = value % p;
  if (r < 0)
    r += p;
  return (int)r;
}

void LCCC_setCoefficientIntegers(void) {
  g_coeff_mode = KH_COEFF_Z;
  g_coeff_modulus = 0;
}

bool LCCC_setCoefficientModPrime(int p) {
  if (!is_prime_int(p))
    return false;

  g_coeff_mode = KH_COEFF_FP;
  g_coeff_modulus = p;
  return true;
}

KHCoefficientMode LCCC_getCoefficientMode(void) {
  return g_coeff_mode;
}

int LCCC_getCoefficientModulus(void) {
  return g_coeff_modulus;
}

int LCCC_coeffNormalize(int value) {
  return normalize_i64((int64_t)value);
}

int LCCC_coeffAdd(int a, int b) {
  return normalize_i64((int64_t)a + (int64_t)b);
}

int LCCC_coeffMultiply(int a, int b) {
  return normalize_i64((int64_t)a * (int64_t)b);
}

int LCCC_coeffNegate(int value) {
  return normalize_i64(-(int64_t)value);
}

bool LCCC_coeffIsZero(int value) {
  return LCCC_coeffNormalize(value) == 0;
}

bool LCCC_coeffIsUnit(int value) {
  int a = LCCC_coeffNormalize(value);
  if (g_coeff_mode == KH_COEFF_Z)
    return a == 1 || a == -1;

  /* F_p is a field because the setter only accepts prime p. */
  return a != 0;
}

bool LCCC_coeffInverse(int value, int *inverse_out) {
  if (inverse_out == NULL)
    return false;

  int a = LCCC_coeffNormalize(value);

  if (g_coeff_mode == KH_COEFF_Z) {
    if (a == 1 || a == -1) {
      *inverse_out = a;
      return true;
    }
    return false;
  }

  if (a == 0)
    return false;

  /* Extended Euclidean algorithm. For prime p and a != 0, gcd(a,p)=1. */
  int64_t old_r = a;
  int64_t r = g_coeff_modulus;
  int64_t old_s = 1;
  int64_t s = 0;

  while (r != 0) {
    int64_t q = old_r / r;

    int64_t tmp = old_r - q * r;
    old_r = r;
    r = tmp;

    tmp = old_s - q * s;
    old_s = s;
    s = tmp;
  }

  if (old_r != 1)
    return false;

  *inverse_out = normalize_i64(old_s);
  return true;
}

/* ================================================================
 *  Internal helpers
 * ================================================================ */

/*
 *  Purpose:    Allocate a new LCCCTerm node.
 */
static LCCCTerm *create_term(CannedCobordism *cc, int coeff) {
  LCCCTerm *t = (LCCCTerm *)malloc(sizeof(LCCCTerm));
  if (t == NULL)
    return NULL;
  t->coeff = LCCC_coeffNormalize(coeff);
  t->cobordism = cc;
  t->next = NULL;
  return t;
}

/*
 *  Purpose:    Check if two cobordisms should be combined as like terms.
 *              We compare by pointer identity first (fast path), then
 *              fall back to structural equality via the CannedCobordism
 *              interface.
 */
static bool cobordism_equal(const CannedCobordism *a,
                            const CannedCobordism *b) {
  if (a == b)
    return true;
  if (a == NULL || b == NULL)
    return false;

  //Structural equality: same source, target, and impl_data
  if (!Cap_equals(a->source, b->source))
    return false;
  if (!Cap_equals(a->target, b->target))
    return false;

  if (a->impl_data == b->impl_data)
    return true;
  if (a->impl_data == NULL || b->impl_data == NULL)
    return false;

  return CannedCobordismImpl_equals((const CannedCobordismImplData *)a->impl_data,
                                    (const CannedCobordismImplData *)b->impl_data);
}

static void lccc_add_term(LCCC *lc, CannedCobordism *cc, int coeff) {
  if (lc == NULL || cc == NULL)
    return;

  coeff = LCCC_coeffNormalize(coeff);
  if (LCCC_coeffIsZero(coeff))
    return;

  LCCCTerm *prev = NULL;
  LCCCTerm *cur = lc->head;

  while (cur != NULL) {
    if (cobordism_equal(cur->cobordism, cc)) {
      cur->coeff = LCCC_coeffAdd(cur->coeff, coeff);

      if (LCCC_coeffIsZero(cur->coeff)) {
        if (prev)
          prev->next = cur->next;
        else
          lc->head = cur->next;

        free(cur);
        lc->count--;
      }

      return;
    }

    prev = cur;
    cur = cur->next;
  }

  LCCCTerm *t = create_term(cc, coeff);
  if (t == NULL)
    return;
  t->next = lc->head;
  lc->head = t;
  lc->count++;
}

/* ================================================================
 *  Construction / Destruction
 * ================================================================ */

LCCC *LCCC_createZero(void) {
  LCCC *lc = (LCCC *)malloc(sizeof(LCCC));
  lc->head = NULL;
  lc->count = 0;
  return lc;
}

LCCC *LCCC_createSingle(CannedCobordism *cc, int coeff) {
  LCCC *lc = LCCC_createZero();
  if (lc == NULL)
    return NULL;

  coeff = LCCC_coeffNormalize(coeff);
  if (cc != NULL && !LCCC_coeffIsZero(coeff)) {
    LCCCTerm *t = create_term(cc, coeff);
    if (t == NULL) {
      free(lc);
      return NULL;
    }
    lc->head = t;
    lc->count = 1;
  }
  return lc;
}

LCCC *LCCC_clone(const LCCC *lc) {
  LCCC *result = LCCC_createZero();
  if (lc == NULL) return result;

  LCCCTerm **tail = &result->head;
  LCCCTerm *cur = lc->head;
  while (cur != NULL) {
    LCCCTerm *t = create_term(cur->cobordism, cur->coeff);
    *tail = t;
    tail = &t->next;
    result->count++;
    cur = cur->next;
  }
  return result;
}

void LCCC_free(LCCC *lc) {
  if (lc == NULL)
    return;
  LCCCTerm *cur = lc->head;
  while (cur != NULL) {
    LCCCTerm *next = cur->next;
    //We do NOT free cur->cobordism -- ownership is external
    free(cur);
    cur = next;
  }
  free(lc);
}

/* ================================================================
 *  Arithmetic (Z or F_p)
 * ================================================================ */

LCCC *LCCC_add(LCCC *a, LCCC *b) {
  LCCC *result = LCCC_clone(a);
  if (b == NULL) return result;

  LCCCTerm *cur = b->head;
  while (cur != NULL) {
    lccc_add_term(result, cur->cobordism, cur->coeff);
    cur = cur->next;
  }
  return result;
}

LCCC *LCCC_compose(LCCC *a, LCCC *b) {
  LCCC *result = LCCC_createZero();
  if (a == NULL || b == NULL || a->count == 0 || b->count == 0) return result;

  LCCCTerm *ai = a->head;
  while (ai != NULL) {
    LCCCTerm *bj = b->head;
    while (bj != NULL) {
      if (ai->cobordism != NULL && bj->cobordism != NULL) {
        CannedCobordism *composed = CannedCobordism_compose(ai->cobordism, bj->cobordism);
        if (composed != NULL) {
          int new_coeff = LCCC_coeffMultiply(ai->coeff, bj->coeff);
          lccc_add_term(result, composed, new_coeff);
        }
      }
      bj = bj->next;
    }
    ai = ai->next;
  }
  return result;
}

LCCC *LCCC_multiply(LCCC *a, RingElement *coeff) {
  if (a == NULL || coeff == NULL || LCCC_coeffIsZero(coeff->value))
    return LCCC_createZero();

  LCCC *result = LCCC_clone(a);
  if (result == NULL)
    return NULL;

  LCCCTerm *cur = result->head;
  while (cur != NULL) {
    cur->coeff = LCCC_coeffMultiply(cur->coeff, coeff->value);
    cur = cur->next;
  }
  return result;
}

/* ================================================================
 * Bar-Natan surface relations
 * ================================================================ */

LCCC *LCCC_reduce(LCCC *a) {
  LCCC *ret = LCCC_createZero();
  if (a == NULL)
    return ret;

  for (LCCCTerm *cur = a->head; cur != NULL; cur = cur->next) {
    if (cur->cobordism == NULL || cur->cobordism->impl_data == NULL)
      continue;

    CannedCobordismImplData *impl =
        (CannedCobordismImplData *)cur->cobordism->impl_data;

    if (!impl->reverse_maps_done) {
      CannedCobordismImpl_reverseMaps(impl);
    }

    int coeff_value = LCCC_coeffNormalize(cur->coeff);
    bool kill = false;

    int nbc = impl->nbc;
    int ncc = impl->ncc;

    int *base_dots = (int *)calloc((size_t)nbc, sizeof(int));
    int *more_work = (int *)malloc((size_t)ncc * sizeof(int));
    if (base_dots == NULL || more_work == NULL) {
      free(base_dots);
      free(more_work);
      continue;
    }

    int nmore = 0;

    /*
     * JavaKh-style reduction:
     *
     * Closed components:
     *   sphere = 0
     *   dotted sphere = 1
     *   torus = 2
     *
     * Components with one boundary:
     *   move dot/genus data to that boundary component.
     *
     * Components with multiple boundaries:
     *   apply neck-cutting. The normal form has each boundary component
     *   as its own connected component, with dots distributed according
     *   to the neck-cutting expansion.
     */
    for (int ci = 0; ci < ncc; ci++) {
      int g = impl->genus[ci];
      int d = impl->dots[ci];
      int bcount = impl->bc_sizes[ci];

      if (g + d > 1) {
        kill = true;
        break;
      }

      if (bcount == 0) {
        if (g == 1 && d == 0) {
          coeff_value = LCCC_coeffMultiply(coeff_value, 2);
        } else if (g == 0 && d == 0) {
          kill = true;
          break;
        } else if (g == 0 && d == 1) {
          //dotted sphere = 1
        } else {
          kill = true;
          break;
        }
      } else if (bcount == 1) {
        int bc = impl->boundaryComponents[ci][0];

        if (bc < 0 || bc >= nbc) {
          kill = true;
          break;
        }

        base_dots[bc] = d + g;

        if (g == 1) {
          coeff_value = LCCC_coeffMultiply(coeff_value, 2);
        }
      } else {
        if (g + d == 1) {
          if (g == 1) {
            coeff_value = LCCC_coeffMultiply(coeff_value, 2);
          }

          for (int j = 0; j < bcount; j++) {
            int bc = impl->boundaryComponents[ci][j];
            if (bc < 0 || bc >= nbc) {
              kill = true;
              break;
            }
            base_dots[bc] = 1;
          }

          if (kill)
            break;
        } else {
          /*
           * Genus 0, dot 0, multiple boundary components.
           * we expand using the neck-cutting relation.
           */
          more_work[nmore++] = ci;
        }
      }
    }

    if (kill || LCCC_coeffIsZero(coeff_value)) {
      free(base_dots);
      free(more_work);
      continue;
    }

    /*
     * Build the neck-cutting dot-pattern expansion.
     *
     * For a component with boundary components b_0,...,b_{m-1},
     * create m terms. In term k, all those boundaries get dots
     * except b_k.
     */
    int pattern_count = 1;
    int **patterns = (int **)malloc(sizeof(int *));
    if (patterns == NULL) {
      free(base_dots);
      free(more_work);
      continue;
    }

    patterns[0] = (int *)malloc((size_t)nbc * sizeof(int));
    if (patterns[0] == NULL) {
      free(patterns);
      free(base_dots);
      free(more_work);
      continue;
    }

    memcpy(patterns[0], base_dots, (size_t)nbc * sizeof(int));

    for (int mw = 0; mw < nmore; mw++) {
      int ci = more_work[mw];
      int bcount = impl->bc_sizes[ci];

      int new_count = pattern_count * bcount;
      int **new_patterns = (int **)malloc((size_t)new_count * sizeof(int *));
      if (new_patterns == NULL) {
        kill = true;
        break;
      }

      int out_idx = 0;

      for (int p = 0; p < pattern_count; p++) {
        for (int choice = 0; choice < bcount; choice++) {
          int *dots = (int *)malloc((size_t)nbc * sizeof(int));
          if (dots == NULL) {
            kill = true;
            break;
          }

          memcpy(dots, patterns[p], (size_t)nbc * sizeof(int));

          for (int j = 0; j < bcount; j++) {
            int bc = impl->boundaryComponents[ci][j];
            if (bc < 0 || bc >= nbc) {
              kill = true;
              free(dots);
              break;
            }
            dots[bc] = 1;
          }

          if (kill)
            break;

          int undotted_bc = impl->boundaryComponents[ci][choice];
          if (undotted_bc < 0 || undotted_bc >= nbc) {
            kill = true;
            free(dots);
            break;
          }

          dots[undotted_bc] = 0;
          new_patterns[out_idx++] = dots;
        }

        if (kill)
          break;
      }

      for (int p = 0; p < pattern_count; p++) {
        free(patterns[p]);
      }
      free(patterns);

      if (kill) {
        for (int p = 0; p < out_idx; p++) {
          free(new_patterns[p]);
        }
        free(new_patterns);
        break;
      }

      patterns = new_patterns;
      pattern_count = new_count;
    }

    if (!kill) {
      for (int p = 0; p < pattern_count; p++) {
        CannedCobordismImplData *new_impl =
            CannedCobordismImpl_create(impl->top, impl->bottom);

        if (new_impl == NULL)
          continue;

        new_impl->hpower = impl->hpower;
        new_impl->ncc = new_impl->nbc;

        for (int bc = 0; bc < new_impl->nbc; bc++) {
          new_impl->connectedComponent[bc] = bc;
        }

        new_impl->dots = (int *)calloc((size_t)new_impl->ncc, sizeof(int));
        new_impl->genus = (int *)calloc((size_t)new_impl->ncc, sizeof(int));

        if (new_impl->dots == NULL || new_impl->genus == NULL) {
          CannedCobordismImpl_free(new_impl);
          continue;
        }

        for (int bc = 0; bc < new_impl->nbc; bc++) {
          new_impl->dots[bc] = patterns[p][bc];
          new_impl->genus[bc] = 0;
        }

        CannedCobordism *new_cc = CannedCobordismImpl_as_CannedCobordism(new_impl);
        lccc_add_term(ret, new_cc, coeff_value);
      }
    }

    for (int p = 0; p < pattern_count; p++) {
      free(patterns[p]);
    }
    free(patterns);
    free(base_dots);
    free(more_work);
  }

  return ret;
}

LCCC *LCCC_invert(LCCC *lc) {
  if (lc == NULL || lc->count != 1 || lc->head == NULL)
    return NULL;

  /*
   * The topological pivot is checked by Komplex before this function is
   * called. Here we only invert its scalar coefficient. Over Z, only +/-1
   * are units. Over F_p, every nonzero coefficient has an inverse.
   */
  int inverse = 0;
  if (!LCCC_coeffInverse(lc->head->coeff, &inverse))
    return NULL;

  return LCCC_createSingle(lc->head->cobordism, inverse);
}

LCCC *LCCC_negate(LCCC *lc) {
  LCCC *result = LCCC_clone(lc);
  if (result == NULL)
    return NULL;

  LCCCTerm *cur = result->head;
  while (cur != NULL) {
    cur->coeff = LCCC_coeffNegate(cur->coeff);
    cur = cur->next;
  }
  return result;
}

LCCC *LCCC_capOffTop(LCCC *lc, Cap *new_top, bool add_dot) {
  if (lc == NULL) return LCCC_createZero();
  LCCC *res = LCCC_createZero();

  /*
   * Accumulate directly into res. The previous implementation created a
   * one-term LCCC and then cloned the entire accumulated list through
   * LCCC_add on every iteration.
   */
  for (LCCCTerm *cur = lc->head; cur != NULL; cur = cur->next) {
    CannedCobordism *capped_cc =
        CannedCobordismImpl_capOffTop(cur->cobordism, new_top, add_dot);
    if (capped_cc != NULL)
      lccc_add_term(res, capped_cc, cur->coeff);
  }

  LCCC *reduced_res = LCCC_reduce(res);
  LCCC_free(res);

  return reduced_res;
}

LCCC *LCCC_cupOnBottom(LCCC *lc, Cap *new_bottom, bool add_dot) {
  if (lc == NULL) return LCCC_createZero();
  LCCC *res = LCCC_createZero();

  for (LCCCTerm *cur = lc->head; cur != NULL; cur = cur->next) {
    CannedCobordism *cupped_cc =
        CannedCobordismImpl_cupOnBottom(cur->cobordism, new_bottom, add_dot);
    if (cupped_cc != NULL)
      lccc_add_term(res, cupped_cc, cur->coeff);
  }

  LCCC *reduced_res = LCCC_reduce(res);
  LCCC_free(res);

  return reduced_res;
}

/* ================================================================
 *  Predicates
 * ================================================================ */

bool LCCC_isZero(LCCC *lc) {
  if (lc == NULL)
    return true;
  return (lc->count == 0);
}

bool LCCC_contains(const LCCC *lc, const CannedCobordism *cc) {
  if (lc == NULL || cc == NULL)
    return false;
  LCCCTerm *cur = lc->head;
  while (cur != NULL) {
    if (cobordism_equal(cur->cobordism, cc))
      return true;
    cur = cur->next;
  }
  return false;
}
