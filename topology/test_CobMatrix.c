/*
 *  test_CobMatrix.c
 *
 *  Purpose:
 *      Standalone test harness for CobMatrix.h / CobMatrix.c.
 *      Provides minimal stub implementations for LCCC, RingElement,
 *      and SmoothingColumn that CobMatrix depends on.
 *
 *  Compile:
 *      gcc -Wall -Wextra -g -o test_CobMatrix.exe \
 *          test_CobMatrix.c CobMatrix.c Cap.c
 *
 *  Run:
 *      ./test_CobMatrix.exe
 */

#include "Cap.h"
#include "CobMatrix.h"


#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 *  Stub implementations for LCCC, RingElement, SmoothingColumn
 *
 *  LCCC is modelled as a simple integer wrapper. This is enough
 *  to exercise the matrix operations (put, add, compose, etc.)
 *  without the real algebra layer.
 * ================================================================ */

struct LCCC {
  int value;
};

struct RingElement {
  int value;
};

LCCC *LCCC_create(int value) {
  LCCC *lc = (LCCC *)malloc(sizeof(LCCC));
  lc->value = value;
  return lc;
}

LCCC *LCCC_add(LCCC *a, LCCC *b) { return LCCC_create(a->value + b->value); }

LCCC *LCCC_compose(LCCC *a, LCCC *b) {
  return LCCC_create(a->value * b->value);
}

LCCC *LCCC_multiply(LCCC *a, RingElement *coeff) {
  return LCCC_create(a->value * coeff->value);
}

LCCC *LCCC_reduce(LCCC *a) { return LCCC_create(a->value); }

bool LCCC_isZero(LCCC *a) { return a->value == 0; }

void LCCC_free(LCCC *a) { free(a); }

SmoothingColumn *SmoothingColumn_clone(SmoothingColumn *col) {
  SmoothingColumn *clone = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  clone->n = col->n;
  clone->numbers = (int *)malloc((size_t)col->n * sizeof(int));
  clone->smoothings = (Cap **)malloc((size_t)col->n * sizeof(Cap *));
  memcpy(clone->numbers, col->numbers, (size_t)col->n * sizeof(int));
  memcpy(clone->smoothings, col->smoothings, (size_t)col->n * sizeof(Cap *));
  return clone;
}

bool SmoothingColumn_equals(SmoothingColumn *a, SmoothingColumn *b) {
  if (a->n != b->n)
    return false;
  for (int i = 0; i < a->n; i++) {
    if (a->numbers[i] != b->numbers[i])
      return false;
  }
  return true;
}

/* ================================================================
 *  Helper: make a simple SmoothingColumn
 * ================================================================ */
static SmoothingColumn *make_smoothing_column(int n) {
  SmoothingColumn *sc = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  sc->n = n;
  sc->numbers = (int *)calloc((size_t)n, sizeof(int));
  sc->smoothings = (Cap **)calloc((size_t)n, sizeof(Cap *));
  for (int i = 0; i < n; i++) {
    sc->numbers[i] = i;
    /* Create tiny caps just for the smoothings pointers */
    int p[] = {1, 0};
    sc->smoothings[i] = Cap_create(2, 0);
    memcpy(sc->smoothings[i]->pairings, p, 2 * sizeof(int));
  }
  return sc;
}

static void free_smoothing_column(SmoothingColumn *sc) {
  for (int i = 0; i < sc->n; i++)
    Cap_free(sc->smoothings[i]);
  free(sc->numbers);
  free(sc->smoothings);
  free(sc);
}

/* ================================================================
 *  Test 1: CobMatrix_create — basic construction (shared mode)
 * ================================================================ */
static void test_create_shared(void) {
  SmoothingColumn *src = make_smoothing_column(3);
  SmoothingColumn *tgt = make_smoothing_column(2);

  CobMatrix *m = CobMatrix_create(src, tgt, true);
  assert(m != NULL);
  assert(m->source == src);
  assert(m->target == tgt);

  /* Empty matrix: isZero should be true */
  assert(CobMatrix_isZero(m) == true);

  CobMatrix_free(m);
  free_smoothing_column(src);
  free_smoothing_column(tgt);
  printf("  [PASS] test_create_shared\n");
}

/* ================================================================
 *  Test 2: CobMatrix_create — cloned mode
 * ================================================================ */
static void test_create_cloned(void) {
  SmoothingColumn *src = make_smoothing_column(3);
  SmoothingColumn *tgt = make_smoothing_column(2);

  CobMatrix *m = CobMatrix_create(src, tgt, false);
  assert(m != NULL);
  /* In non-shared mode, source/target are clones */
  assert(m->source != src);
  assert(m->target != tgt);
  assert(m->source->n == src->n);
  assert(m->target->n == tgt->n);

  CobMatrix_free(m);
  free_smoothing_column(src);
  free_smoothing_column(tgt);
  printf("  [PASS] test_create_cloned\n");
}

/* ================================================================
 *  Test 3: putEntry and isZero
 * ================================================================ */
static void test_putEntry(void) {
  SmoothingColumn *src = make_smoothing_column(3);
  SmoothingColumn *tgt = make_smoothing_column(2);

  CobMatrix *m = CobMatrix_create(src, tgt, true);
  assert(CobMatrix_isZero(m) == true);

  LCCC *val = LCCC_create(5);
  CobMatrix_putEntry(m, 0, 1, val);
  assert(CobMatrix_isZero(m) == false);

  CobMatrix_free(m);
  free_smoothing_column(src);
  free_smoothing_column(tgt);
  printf("  [PASS] test_putEntry\n");
}

/* ================================================================
 *  Test 4: putEntry — NULL / zero value is ignored
 * ================================================================ */
static void test_putEntry_null(void) {
  SmoothingColumn *src = make_smoothing_column(2);
  SmoothingColumn *tgt = make_smoothing_column(2);

  CobMatrix *m = CobMatrix_create(src, tgt, true);

  CobMatrix_putEntry(m, 0, 0, NULL);
  assert(CobMatrix_isZero(m) == true);

  LCCC *zero = LCCC_create(0);
  CobMatrix_putEntry(m, 0, 0, zero);
  assert(CobMatrix_isZero(m) == true);

  LCCC_free(zero);
  CobMatrix_free(m);
  free_smoothing_column(src);
  free_smoothing_column(tgt);
  printf("  [PASS] test_putEntry_null\n");
}

/* ================================================================
 *  Test 5: unpackRow
 * ================================================================ */
static void test_unpackRow(void) {
  SmoothingColumn *src = make_smoothing_column(3);
  SmoothingColumn *tgt = make_smoothing_column(2);

  CobMatrix *m = CobMatrix_create(src, tgt, true);

  LCCC *v0 = LCCC_create(10);
  LCCC *v2 = LCCC_create(20);
  CobMatrix_putEntry(m, 1, 0, v0);
  CobMatrix_putEntry(m, 1, 2, v2);

  LCCC **row = CobMatrix_unpackRow(m, 1);
  assert(row[0] != NULL && row[0]->value == 10);
  assert(row[1] == NULL);
  assert(row[2] != NULL && row[2]->value == 20);

  free(row);
  CobMatrix_free(m);
  free_smoothing_column(src);
  free_smoothing_column(tgt);
  printf("  [PASS] test_unpackRow\n");
}

/* ================================================================
 *  Test 6: addEntry — new entry
 * ================================================================ */
static void test_addEntry_new(void) {
  SmoothingColumn *src = make_smoothing_column(2);
  SmoothingColumn *tgt = make_smoothing_column(2);

  CobMatrix *m = CobMatrix_create(src, tgt, true);

  LCCC *val = LCCC_create(7);
  CobMatrix_addEntry(m, 0, 0, val);
  assert(CobMatrix_isZero(m) == false);

  LCCC **row = CobMatrix_unpackRow(m, 0);
  assert(row[0] != NULL && row[0]->value == 7);

  free(row);
  CobMatrix_free(m);
  free_smoothing_column(src);
  free_smoothing_column(tgt);
  printf("  [PASS] test_addEntry_new\n");
}

/* ================================================================
 *  Test 7: addEntry — adding to existing entry
 * ================================================================ */
static void test_addEntry_accumulate(void) {
  SmoothingColumn *src = make_smoothing_column(2);
  SmoothingColumn *tgt = make_smoothing_column(2);

  CobMatrix *m = CobMatrix_create(src, tgt, true);

  LCCC *v1 = LCCC_create(3);
  LCCC *v2 = LCCC_create(4);
  CobMatrix_addEntry(m, 0, 0, v1);
  CobMatrix_addEntry(m, 0, 0, v2);

  LCCC **row = CobMatrix_unpackRow(m, 0);
  assert(row[0] != NULL && row[0]->value == 7); /* 3 + 4 */

  free(row);
  CobMatrix_free(m);
  free_smoothing_column(src);
  free_smoothing_column(tgt);
  printf("  [PASS] test_addEntry_accumulate\n");
}

/* ================================================================
 *  Test 8: addEntry — cancellation (sum to zero removes entry)
 * ================================================================ */
static void test_addEntry_cancel(void) {
  SmoothingColumn *src = make_smoothing_column(2);
  SmoothingColumn *tgt = make_smoothing_column(2);

  CobMatrix *m = CobMatrix_create(src, tgt, true);

  LCCC *v1 = LCCC_create(5);
  LCCC *v2 = LCCC_create(-5);
  CobMatrix_addEntry(m, 0, 0, v1);
  CobMatrix_addEntry(m, 0, 0, v2);

  assert(CobMatrix_isZero(m) == true);

  CobMatrix_free(m);
  free_smoothing_column(src);
  free_smoothing_column(tgt);
  printf("  [PASS] test_addEntry_cancel\n");
}

/* ================================================================
 *  Test 9: compose — identity-like (1 on diagonal)
 *
 *  m1 = [[1,0],[0,1]], m2 = [[a,b],[c,d]]
 *  m1 * m2 should give [[a,b],[c,d]]
 * ================================================================ */
static void test_compose_identity(void) {
  SmoothingColumn *sc = make_smoothing_column(2);

  /* m1 = identity matrix */
  CobMatrix *m1 = CobMatrix_create(sc, sc, true);
  CobMatrix_putEntry(m1, 0, 0, LCCC_create(1));
  CobMatrix_putEntry(m1, 1, 1, LCCC_create(1));

  /* m2 = [[2,3],[4,5]] */
  CobMatrix *m2 = CobMatrix_create(sc, sc, true);
  CobMatrix_putEntry(m2, 0, 0, LCCC_create(2));
  CobMatrix_putEntry(m2, 0, 1, LCCC_create(3));
  CobMatrix_putEntry(m2, 1, 0, LCCC_create(4));
  CobMatrix_putEntry(m2, 1, 1, LCCC_create(5));

  CobMatrix *result = CobMatrix_compose(m1, m2);
  assert(result != NULL);
  assert(CobMatrix_isZero(result) == false);

  /* Check result[0][0] = 1*2 + 0*4 = 2 */
  LCCC **row0 = CobMatrix_unpackRow(result, 0);
  assert(row0[0] != NULL && row0[0]->value == 2);
  assert(row0[1] != NULL && row0[1]->value == 3);
  free(row0);

  /* Check result[1][0] = 0*2 + 1*4 = 4 */
  LCCC **row1 = CobMatrix_unpackRow(result, 1);
  assert(row1[0] != NULL && row1[0]->value == 4);
  assert(row1[1] != NULL && row1[1]->value == 5);
  free(row1);

  CobMatrix_free(result);
  CobMatrix_free(m1);
  CobMatrix_free(m2);
  free_smoothing_column(sc);
  printf("  [PASS] test_compose_identity\n");
}

/* ================================================================
 *  Test 10: compose — non-trivial matrix multiplication
 *
 *  m1 = [[1,2],[3,4]], m2 = [[5,6],[7,8]]
 *  result = m1 * m2 = [[1*5+2*7, 1*6+2*8], [3*5+4*7, 3*6+4*8]]
 *                    = [[19, 22], [43, 50]]
 * ================================================================ */
static void test_compose_nontrivial(void) {
  SmoothingColumn *sc = make_smoothing_column(2);

  CobMatrix *m1 = CobMatrix_create(sc, sc, true);
  CobMatrix_putEntry(m1, 0, 0, LCCC_create(1));
  CobMatrix_putEntry(m1, 0, 1, LCCC_create(2));
  CobMatrix_putEntry(m1, 1, 0, LCCC_create(3));
  CobMatrix_putEntry(m1, 1, 1, LCCC_create(4));

  CobMatrix *m2 = CobMatrix_create(sc, sc, true);
  CobMatrix_putEntry(m2, 0, 0, LCCC_create(5));
  CobMatrix_putEntry(m2, 0, 1, LCCC_create(6));
  CobMatrix_putEntry(m2, 1, 0, LCCC_create(7));
  CobMatrix_putEntry(m2, 1, 1, LCCC_create(8));

  CobMatrix *result = CobMatrix_compose(m1, m2);

  LCCC **row0 = CobMatrix_unpackRow(result, 0);
  assert(row0[0] != NULL && row0[0]->value == 19);
  assert(row0[1] != NULL && row0[1]->value == 22);
  free(row0);

  LCCC **row1 = CobMatrix_unpackRow(result, 1);
  assert(row1[0] != NULL && row1[0]->value == 43);
  assert(row1[1] != NULL && row1[1]->value == 50);
  free(row1);

  CobMatrix_free(result);
  CobMatrix_free(m1);
  CobMatrix_free(m2);
  free_smoothing_column(sc);
  printf("  [PASS] test_compose_nontrivial\n");
}

/* ================================================================
 *  Test 11: multiply — scalar multiplication
 * ================================================================ */
static void test_multiply(void) {
  SmoothingColumn *sc = make_smoothing_column(2);

  CobMatrix *m = CobMatrix_create(sc, sc, true);
  CobMatrix_putEntry(m, 0, 0, LCCC_create(3));
  CobMatrix_putEntry(m, 1, 1, LCCC_create(5));

  RingElement coeff = {2};
  CobMatrix_multiply(m, &coeff);

  LCCC **row0 = CobMatrix_unpackRow(m, 0);
  assert(row0[0] != NULL && row0[0]->value == 6); /* 3*2 */
  free(row0);

  LCCC **row1 = CobMatrix_unpackRow(m, 1);
  assert(row1[1] != NULL && row1[1]->value == 10); /* 5*2 */
  free(row1);

  CobMatrix_free(m);
  free_smoothing_column(sc);
  printf("  [PASS] test_multiply\n");
}

/* ================================================================
 *  Test 12: add — matrix addition
 * ================================================================ */
static void test_add(void) {
  SmoothingColumn *sc = make_smoothing_column(2);

  CobMatrix *dest = CobMatrix_create(sc, sc, true);
  CobMatrix_putEntry(dest, 0, 0, LCCC_create(1));
  CobMatrix_putEntry(dest, 1, 1, LCCC_create(2));

  CobMatrix *src = CobMatrix_create(sc, sc, true);
  CobMatrix_putEntry(src, 0, 0, LCCC_create(10));
  CobMatrix_putEntry(src, 0, 1, LCCC_create(20));

  CobMatrix_add(dest, src);

  LCCC **row0 = CobMatrix_unpackRow(dest, 0);
  assert(row0[0] != NULL && row0[0]->value == 11); /* 1+10 */
  assert(row0[1] != NULL && row0[1]->value == 20); /* 0+20 */
  free(row0);

  CobMatrix_free(dest);
  CobMatrix_free(src);
  free_smoothing_column(sc);
  printf("  [PASS] test_add\n");
}

/* ================================================================
 *  Test 13: reduce — entries are reduced (our stub is identity)
 * ================================================================ */
static void test_reduce(void) {
  SmoothingColumn *sc = make_smoothing_column(2);

  CobMatrix *m = CobMatrix_create(sc, sc, true);
  CobMatrix_putEntry(m, 0, 0, LCCC_create(7));

  CobMatrix_reduce(m);
  assert(CobMatrix_isZero(m) == false);

  LCCC **row0 = CobMatrix_unpackRow(m, 0);
  assert(row0[0] != NULL && row0[0]->value == 7);
  free(row0);

  CobMatrix_free(m);
  free_smoothing_column(sc);
  printf("  [PASS] test_reduce\n");
}

/* ================================================================
 *  Test 14: isZero — empty matrix
 * ================================================================ */
static void test_isZero_empty(void) {
  SmoothingColumn *sc = make_smoothing_column(3);
  CobMatrix *m = CobMatrix_create(sc, sc, true);

  assert(CobMatrix_isZero(m) == true);

  CobMatrix_free(m);
  free_smoothing_column(sc);
  printf("  [PASS] test_isZero_empty\n");
}

/* ================================================================
 *  Test 15: isZero — non-empty matrix
 * ================================================================ */
static void test_isZero_nonempty(void) {
  SmoothingColumn *sc = make_smoothing_column(2);
  CobMatrix *m = CobMatrix_create(sc, sc, true);

  CobMatrix_putEntry(m, 0, 0, LCCC_create(1));
  assert(CobMatrix_isZero(m) == false);

  CobMatrix_free(m);
  free_smoothing_column(sc);
  printf("  [PASS] test_isZero_nonempty\n");
}

/* ================================================================
 *  Test 16: extractColumn
 * ================================================================ */
static void test_extractColumn(void) {
  SmoothingColumn *src = make_smoothing_column(3);
  SmoothingColumn *tgt = make_smoothing_column(2);

  CobMatrix *m = CobMatrix_create(src, tgt, true);
  CobMatrix_putEntry(m, 0, 0, LCCC_create(1));
  CobMatrix_putEntry(m, 0, 1, LCCC_create(2));
  CobMatrix_putEntry(m, 0, 2, LCCC_create(3));
  CobMatrix_putEntry(m, 1, 1, LCCC_create(4));

  /* Extract column 1 */
  CobMatrix *col = CobMatrix_extractColumn(m, 1);
  assert(col != NULL);
  assert(col->source->n == 1);

  /* The extracted column should have the value from row 0, col 1 */
  LCCC **col_row0 = CobMatrix_unpackRow(col, 0);
  assert(col_row0[0] != NULL && col_row0[0]->value == 2);
  free(col_row0);

  LCCC **col_row1 = CobMatrix_unpackRow(col, 1);
  assert(col_row1[0] != NULL && col_row1[0]->value == 4);
  free(col_row1);

  /* Original matrix should now have source->n == 2 */
  assert(m->source->n == 2);

  CobMatrix_free(col);
  CobMatrix_free(m);
  free_smoothing_column(src);
  free_smoothing_column(tgt);
  printf("  [PASS] test_extractColumn\n");
}

/* ================================================================
 *  Test 17: extractRow
 * ================================================================ */
static void test_extractRow(void) {
  SmoothingColumn *src = make_smoothing_column(2);
  SmoothingColumn *tgt = make_smoothing_column(3);

  CobMatrix *m = CobMatrix_create(src, tgt, true);
  CobMatrix_putEntry(m, 0, 0, LCCC_create(10));
  CobMatrix_putEntry(m, 1, 0, LCCC_create(20));
  CobMatrix_putEntry(m, 1, 1, LCCC_create(30));
  CobMatrix_putEntry(m, 2, 1, LCCC_create(40));

  /* Extract row 1 */
  CobMatrix *row = CobMatrix_extractRow(m, 1);
  assert(row != NULL);
  assert(row->target->n == 1);

  /* The extracted row should have values from original row 1 */
  LCCC **row_row0 = CobMatrix_unpackRow(row, 0);
  assert(row_row0[0] != NULL && row_row0[0]->value == 20);
  assert(row_row0[1] != NULL && row_row0[1]->value == 30);
  free(row_row0);

  /* Original matrix should now have target->n == 2 */
  assert(m->target->n == 2);

  CobMatrix_free(row);
  CobMatrix_free(m);
  free_smoothing_column(src);
  free_smoothing_column(tgt);
  printf("  [PASS] test_extractRow\n");
}

/* ================================================================
 *  Test 18: free — NULL is safe
 * ================================================================ */
static void test_free_null(void) {
  CobMatrix_free(NULL);
  printf("  [PASS] test_free_null\n");
}

/* ================================================================
 *  main
 * ================================================================ */
int main(void) {
  printf("Running CobMatrix tests...\n\n");

  test_create_shared();
  test_create_cloned();
  test_putEntry();
  test_putEntry_null();
  test_unpackRow();
  test_addEntry_new();
  test_addEntry_accumulate();
  test_addEntry_cancel();
  test_compose_identity();
  test_compose_nontrivial();
  test_multiply();
  test_add();
  test_reduce();
  test_isZero_empty();
  test_isZero_nonempty();
  test_extractColumn();
  test_extractRow();
  test_free_null();

  printf("\nAll CobMatrix tests passed!\n");
  return 0;
}
