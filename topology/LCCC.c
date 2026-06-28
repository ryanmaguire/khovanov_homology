/*
 *  LCCC.c
 *
 *  Linear Combination of Canned Cobordisms with integer coefficients.
 *
 *  Purpose:
 *      Implements the coefficient type used in the cobordism-valued
 *      Khovanov complex. Like cobordism terms are combined by adding
 *      integer coefficients, and reduction applies the Bar-Natan
 *      surface relations.
 *
 *  Compile (standalone test):
 *      gcc -Wall -Wextra -g -c LCCC.c -I.
 */

#include "LCCC.h"
#include "CannedCobordismImpl.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 *  Internal helpers
 * ================================================================ */

/*
 *  Purpose:    Allocate a new LCCCTerm node.
 */
static LCCCTerm *create_term(CannedCobordism *cc, int coeff) {
  LCCCTerm *t = (LCCCTerm *)malloc(sizeof(LCCCTerm));
  t->coeff = coeff;
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
  if (lc == NULL || cc == NULL || coeff == 0) return;

  LCCCTerm *prev = NULL;
  LCCCTerm *cur = lc->head;

  while (cur != NULL) {
    if (cobordism_equal(cur->cobordism, cc)) {
      long merged = (long)cur->coeff + (long)coeff;

      if (merged > INT_MAX || merged < INT_MIN) {
        abort();  //Coefficient overflow; replace int with int64_t if needed.
      }

      cur->coeff = (int)merged;

      if (cur->coeff == 0) {
        if (prev) prev->next = cur->next;
        else lc->head = cur->next;

        free(cur);
        lc->count--;
      }

      return;
    }

    prev = cur;
    cur = cur->next;
  }

  LCCCTerm *t = create_term(cc, coeff);
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
  if (cc != NULL && coeff != 0) {
    LCCCTerm *t = create_term(cc, coeff);
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
 *  Arithmetic (Z)
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
          int new_coeff = ai->coeff * bj->coeff;
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
  if (a == NULL || coeff == NULL || coeff->value == 0) return LCCC_createZero();
  
  LCCC *result = LCCC_clone(a);
  LCCCTerm *cur = result->head;
  while (cur != NULL) {
      cur->coeff *= coeff->value;
      cur = cur->next;
  }
  return result;
}

/* ================================================================
 * The Rulebook (Bar-Natan Relations)
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

    long coeff_long = (long)cur->coeff;
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
          coeff_long *= 2;
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
          coeff_long *= 2;
        }
      } else {
        if (g + d == 1) {
          if (g == 1) {
            coeff_long *= 2;
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
           * This is the real neck-cutting work.
           */
          more_work[nmore++] = ci;
        }
      }
    }

    if (kill || coeff_long == 0) {
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

        if (coeff_long > INT_MAX || coeff_long < INT_MIN) {
          CannedCobordismImpl_free(new_impl);
          continue;
        }

        CannedCobordism *new_cc = CannedCobordismImpl_as_CannedCobordism(new_impl);
        lccc_add_term(ret, new_cc, (int)coeff_long);
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
   * Only unit pivots with coefficient +/-1 are invertible over Z.
   * The current Gaussian elimination path uses identity-like cobordism
   * pivots, whose inverse is represented by the same underlying term.
   */
  if (lc->head->coeff != 1 && lc->head->coeff != -1)
    return NULL;

  return LCCC_clone(lc);
}

LCCC *LCCC_negate(LCCC *lc) {
  LCCC *result = LCCC_clone(lc);
  LCCCTerm *cur = result->head;
  while (cur != NULL) {
      cur->coeff *= -1;
      cur = cur->next;
  }
  return result;
}

LCCC *LCCC_capOffTop(LCCC *lc, Cap *new_top, bool add_dot) {
    if (lc == NULL) return LCCC_createZero();
    LCCC *res = LCCC_createZero();

    LCCCTerm *cur = lc->head;
    while (cur != NULL) {
        CannedCobordism *capped_cc = CannedCobordismImpl_capOffTop(cur->cobordism, new_top, add_dot);
        if (capped_cc != NULL) {
            LCCC *single = LCCC_createSingle(capped_cc, cur->coeff);
            LCCC *sum = LCCC_add(res, single);
            LCCC_free(res); 
            LCCC_free(single);
            res = sum;
        }
        cur = cur->next;
    }
    
    LCCC *reduced_res = LCCC_reduce(res);
    LCCC_free(res);
    
    return reduced_res;
}

LCCC *LCCC_cupOnBottom(LCCC *lc, Cap *new_bottom, bool add_dot) {
    if (lc == NULL) return LCCC_createZero();
    LCCC *res = LCCC_createZero();

    LCCCTerm *cur = lc->head;
    while (cur != NULL) {
        CannedCobordism *cupped_cc = CannedCobordismImpl_cupOnBottom(cur->cobordism, new_bottom, add_dot);
        if (cupped_cc != NULL) {
            LCCC *single = LCCC_createSingle(cupped_cc, cur->coeff);
            LCCC *sum = LCCC_add(res, single);
            LCCC_free(res); 
            LCCC_free(single);
            res = sum;
        }
        cur = cur->next;
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
