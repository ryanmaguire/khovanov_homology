/*
 *  test_CannedCobordismImpl.c
 *
 *  Purpose:
 *      Standalone test harness for CannedCobordismImpl.h / .c.
 *      Uses assert() — no external test framework required.
 *
 *  Compile:
 *      gcc -Wall -Wextra -g -o test_CannedCobordismImpl.exe \
 *          test_CannedCobordismImpl.c CannedCobordismImpl.c Cap.c
 *
 *  Run:
 *      .\test_CannedCobordismImpl.exe
 */

#include "CannedCobordism.h"
#include "CannedCobordismImpl.h"
#include "Cap.h"

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 *  Helper: build a Cap with given pairings (heap-allocated)
 * ================================================================ */
static Cap *make_cap(int n, int ncycles, const int *pairings) {
  Cap *c = Cap_create(n, ncycles);
  assert(c != NULL);
  memcpy(c->pairings, pairings, (size_t)n * sizeof(int));
  return c;
}

/* ================================================================
 *  Test 1: create — simple 4-edge, same top & bottom
 *
 *  Pairings {1,0,3,2}: edges 0↔1, 2↔3. Two arcs.
 *  With top == bottom and these pairings, tracing gives:
 *    start at 0 → top pairs 0→1 → bottom pairs 1→0 → done (BC 0)
 *    start at 2 → top pairs 2→3 → bottom pairs 3→2 → done (BC 1)
 *  So offtop = 2 mixed BCs, nbc = 2 (no cycles).
 * ================================================================ */
static void test_create_simple(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);

  CannedCobordismImplData *d = CannedCobordismImpl_create(c, c);
  assert(d != NULL);
  assert(d->n == 4);
  assert(d->offtop == 2); /* 2 mixed boundary components */
  assert(d->offbot == 2); /* no top cycles */
  assert(d->nbc == 2);    /* no bottom cycles either */

  /* component[0] == component[1], component[2] == component[3] */
  assert(d->component[0] == d->component[1]);
  assert(d->component[2] == d->component[3]);
  assert(d->component[0] != d->component[2]);

  CannedCobordismImpl_free(d);
  Cap_free(c);
  printf("  [PASS] test_create_simple\n");
}

/* ================================================================
 *  Test 2: create — different top & bottom pairings
 *
 *  top:    {1,0,3,2}   edges 0↔1, 2↔3
 *  bottom: {3,2,1,0}   edges 0↔3, 1↔2
 *
 *  Tracing from edge 0:
 *    0 → top[0]=1 → bottom[1]=2 → top[2]=3 → bottom[3]=0 → done
 *  All 4 edges in a single boundary component.
 *  offtop = 1, nbc = 1.
 * ================================================================ */
static void test_create_different_caps(void) {
  int pt[] = {1, 0, 3, 2};
  int pb[] = {3, 2, 1, 0};
  Cap *top = make_cap(4, 0, pt);
  Cap *bot = make_cap(4, 0, pb);

  CannedCobordismImplData *d = CannedCobordismImpl_create(top, bot);
  assert(d != NULL);
  assert(d->n == 4);
  assert(d->offtop == 1); /* single mixed BC */
  assert(d->nbc == 1);

  /* All edges in the same component */
  for (int i = 1; i < 4; i++)
    assert(d->component[i] == d->component[0]);

  CannedCobordismImpl_free(d);
  Cap_free(top);
  Cap_free(bot);
  printf("  [PASS] test_create_different_caps\n");
}

/* ================================================================
 *  Test 3: isomorphism factory
 *
 *  isomorphism(cap) should yield a cobordism where:
 *    - top == bottom == cap
 *    - ncc == nbc
 *    - connectedComponent[i] == i
 *    - dots[i] == 0, genus[i] == 0
 *    - isIsomorphism returns true
 *    - check returns true
 * ================================================================ */
static void test_isomorphism_factory(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);

  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  assert(iso != NULL);

  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;
  assert(d->ncc == d->nbc);
  assert(Cap_equals(d->top, d->bottom));

  for (int i = 0; i < d->nbc; i++) {
    assert(d->connectedComponent[i] == i);
  }
  for (int i = 0; i < d->ncc; i++) {
    assert(d->dots[i] == 0);
    assert(d->genus[i] == 0);
  }

  assert(CannedCobordismImpl_isIsomorphism(d) == true);
  assert(CannedCobordismImpl_check(d) == true);

  CannedCobordism_free(iso);
  Cap_free(c);
  printf("  [PASS] test_isomorphism_factory\n");
}

/* ================================================================
 *  Test 4: isIsomorphism — false when dots are nonzero
 * ================================================================ */
static void test_isIsomorphism_false(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);

  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;

  /* Mutate: set a dot → no longer an isomorphism */
  d->dots[0] = 1;
  assert(CannedCobordismImpl_isIsomorphism(d) == false);

  /* Restore and mutate genus instead */
  d->dots[0] = 0;
  d->genus[0] = 1;
  assert(CannedCobordismImpl_isIsomorphism(d) == false);

  /* Restore and mutate hpower */
  d->genus[0] = 0;
  d->hpower = 1;
  assert(CannedCobordismImpl_isIsomorphism(d) == false);

  CannedCobordism_free(iso);
  Cap_free(c);
  printf("  [PASS] test_isIsomorphism_false\n");
}

/* ================================================================
 *  Test 5: check — valid cobordism
 * ================================================================ */
static void test_check_valid(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);
  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;

  assert(CannedCobordismImpl_check(d) == true);

  CannedCobordism_free(iso);
  Cap_free(c);
  printf("  [PASS] test_check_valid\n");
}

/* ================================================================
 *  Test 6: check — invalid: ncc = 0 with nbc > 0
 * ================================================================ */
static void test_check_invalid(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);
  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;

  int8_t saved_ncc = d->ncc;
  d->ncc = 0;
  assert(CannedCobordismImpl_check(d) == false);
  d->ncc = saved_ncc; /* restore before freeing */

  CannedCobordism_free(iso);
  Cap_free(c);
  printf("  [PASS] test_check_invalid\n");
}

/* ================================================================
 *  Test 7: reverseMaps — verify contents
 *
 *  For an isomorphism with pairings {1,0,3,2}:
 *  nbc == ncc == 2, offtop == 2.
 *  BC 0 has edges {0,1}, BC 1 has edges {2,3}.
 *  Each BC is its own CC, so boundaryComponents[i] = {i}.
 * ================================================================ */
static void test_reverseMaps(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);
  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;

  assert(d->reverse_maps_done == false);
  CannedCobordismImpl_reverseMaps(d);
  assert(d->reverse_maps_done == true);

  /* Each CC has exactly 1 boundary component (itself) */
  for (int i = 0; i < d->ncc; i++) {
    assert(d->bc_sizes[i] == 1);
    assert(d->boundaryComponents[i][0] == i);
  }

  /* Each mixed BC has 2 edges */
  for (int i = 0; i < d->offtop; i++) {
    assert(d->edge_sizes[i] == 2);
  }

  CannedCobordism_free(iso);
  Cap_free(c);
  printf("  [PASS] test_reverseMaps\n");
}

/* ================================================================
 *  Test 8: reverseMaps — idempotent (call twice, no crash)
 * ================================================================ */
static void test_reverseMaps_idempotent(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);
  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;

  CannedCobordismImpl_reverseMaps(d);
  int bc0 = d->bc_sizes[0];
  CannedCobordismImpl_reverseMaps(d);
  assert(d->bc_sizes[0] == bc0);
  assert(d->reverse_maps_done == true);

  CannedCobordism_free(iso);
  Cap_free(c);
  printf("  [PASS] test_reverseMaps_idempotent\n");
}

/* ================================================================
 *  Test 9: equals — self
 * ================================================================ */
static void test_equals_self(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);
  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;

  assert(CannedCobordismImpl_equals(d, d) == true);

  CannedCobordism_free(iso);
  Cap_free(c);
  printf("  [PASS] test_equals_self\n");
}

/* ================================================================
 *  Test 10: equals — two identical isomorphisms
 * ================================================================ */
static void test_equals_cloned(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);

  CannedCobordism *iso1 = CannedCobordismImpl_isomorphism(c);
  CannedCobordism *iso2 = CannedCobordismImpl_isomorphism(c);

  CannedCobordismImplData *d1 = (CannedCobordismImplData *)iso1->impl_data;
  CannedCobordismImplData *d2 = (CannedCobordismImplData *)iso2->impl_data;

  assert(CannedCobordismImpl_equals(d1, d2) == true);

  CannedCobordism_free(iso1);
  CannedCobordism_free(iso2);
  Cap_free(c);
  printf("  [PASS] test_equals_cloned\n");
}

/* ================================================================
 *  Test 11: equals — different caps → not equal
 * ================================================================ */
static void test_equals_different(void) {
  int p1[] = {1, 0, 3, 2};
  int p2[] = {3, 2, 1, 0};
  Cap *c1 = make_cap(4, 0, p1);
  Cap *c2 = make_cap(4, 0, p2);

  CannedCobordism *iso1 = CannedCobordismImpl_isomorphism(c1);
  CannedCobordism *iso2 = CannedCobordismImpl_isomorphism(c2);

  CannedCobordismImplData *d1 = (CannedCobordismImplData *)iso1->impl_data;
  CannedCobordismImplData *d2 = (CannedCobordismImplData *)iso2->impl_data;

  assert(CannedCobordismImpl_equals(d1, d2) == false);

  CannedCobordism_free(iso1);
  CannedCobordism_free(iso2);
  Cap_free(c1);
  Cap_free(c2);
  printf("  [PASS] test_equals_different\n");
}

/* ================================================================
 *  Test 12: hashCode — same impl → same hash
 * ================================================================ */
static void test_hashCode_same(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);

  CannedCobordism *iso1 = CannedCobordismImpl_isomorphism(c);
  CannedCobordism *iso2 = CannedCobordismImpl_isomorphism(c);

  CannedCobordismImplData *d1 = (CannedCobordismImplData *)iso1->impl_data;
  CannedCobordismImplData *d2 = (CannedCobordismImplData *)iso2->impl_data;

  assert(CannedCobordismImpl_hashCode(d1) == CannedCobordismImpl_hashCode(d2));

  CannedCobordism_free(iso1);
  CannedCobordism_free(iso2);
  Cap_free(c);
  printf("  [PASS] test_hashCode_same\n");
}

/* ================================================================
 *  Test 13: hashCode — caching (second call returns same value)
 * ================================================================ */
static void test_hashCode_caching(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);
  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;

  int h1 = CannedCobordismImpl_hashCode(d);
  int h2 = CannedCobordismImpl_hashCode(d);
  assert(h1 == h2);
  assert(d->hashcode == h1);

  CannedCobordism_free(iso);
  Cap_free(c);
  printf("  [PASS] test_hashCode_caching\n");
}

/* ================================================================
 *  Test 14: compareTo — self → 0
 * ================================================================ */
static void test_compareTo_self(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);
  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;

  assert(CannedCobordismImpl_compareTo(d, d) == 0);

  CannedCobordism_free(iso);
  Cap_free(c);
  printf("  [PASS] test_compareTo_self\n");
}

/* ================================================================
 *  Test 15: compareTo — ordering: a < b ↔ b > a
 * ================================================================ */
static void test_compareTo_ordering(void) {
  int p1[] = {1, 0, 3, 2};
  int p2[] = {3, 2, 1, 0};
  Cap *c1 = make_cap(4, 0, p1);
  Cap *c2 = make_cap(4, 0, p2);

  CannedCobordism *iso1 = CannedCobordismImpl_isomorphism(c1);
  CannedCobordism *iso2 = CannedCobordismImpl_isomorphism(c2);

  CannedCobordismImplData *d1 = (CannedCobordismImplData *)iso1->impl_data;
  CannedCobordismImplData *d2 = (CannedCobordismImplData *)iso2->impl_data;

  int cmp1 = CannedCobordismImpl_compareTo(d1, d2);
  int cmp2 = CannedCobordismImpl_compareTo(d2, d1);

  /* Opposite signs (or both zero, but these should differ) */
  assert((cmp1 < 0 && cmp2 > 0) || (cmp1 > 0 && cmp2 < 0) ||
         (cmp1 == 0 && cmp2 == 0));
  /* They're different cobordisms, so should not be zero */
  assert(cmp1 != 0);

  CannedCobordism_free(iso1);
  CannedCobordism_free(iso2);
  Cap_free(c1);
  Cap_free(c2);
  printf("  [PASS] test_compareTo_ordering\n");
}

/* ================================================================
 *  Test 16: compose_vertical — identity ∘ identity == identity
 *
 *  Stacking an isomorphism on itself should yield an isomorphism.
 * ================================================================ */
static void test_compose_vertical_identity(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);

  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;

  CannedCobordism *composed = CannedCobordismImpl_compose_vertical(d, d);
  assert(composed != NULL);

  CannedCobordismImplData *rd = (CannedCobordismImplData *)composed->impl_data;
  assert(CannedCobordismImpl_check(rd) == true);
  assert(CannedCobordismImpl_isIsomorphism(rd) == true);

  /* top and bottom of composed should match the original cap */
  assert(Cap_equals(rd->top, c));
  assert(Cap_equals(rd->bottom, c));

  CannedCobordism_free(composed);
  CannedCobordism_free(iso);
  Cap_free(c);
  printf("  [PASS] test_compose_vertical_identity\n");
}

/* ================================================================
 *  Test 17: compose_vertical — non-trivial cobordisms
 *
 *  Build two cobordisms with different top/bottom that share a
 *  common interface, compose them, and verify consistency.
 *
 *  cob1: top={1,0,3,2}  bottom={3,2,1,0}   (merge)
 *  cob2: top={3,2,1,0}  bottom={1,0,3,2}   (split)
 *  cob1.top == cob2.bottom, so we can compose cob2 on top of cob1.
 * ================================================================ */
static void test_compose_vertical_nontrivial(void) {
  int pa[] = {1, 0, 3, 2};
  int pb[] = {3, 2, 1, 0};
  Cap *capA = make_cap(4, 0, pa);
  Cap *capB = make_cap(4, 0, pb);

  /* cob1:  top=capA, bottom=capB */
  CannedCobordismImplData *d1 = CannedCobordismImpl_create(capA, capB);
  /* Set up as single CC: all boundary components in CC 0 */
  d1->ncc = 1;
  for (int i = 0; i < d1->nbc; i++)
    d1->connectedComponent[i] = 0;
  d1->dots = (int8_t *)calloc((size_t)d1->ncc, sizeof(int8_t));
  d1->genus = (int8_t *)calloc((size_t)d1->ncc, sizeof(int8_t));
  d1->dots[0] = 1;
  assert(CannedCobordismImpl_check(d1));

  /* cob2:  top=capB, bottom=capA  (so cob2.bottom == cob1.top) */
  CannedCobordismImplData *d2 = CannedCobordismImpl_create(capB, capA);
  d2->ncc = 1;
  for (int i = 0; i < d2->nbc; i++)
    d2->connectedComponent[i] = 0;
  d2->dots = (int8_t *)calloc((size_t)d2->ncc, sizeof(int8_t));
  d2->genus = (int8_t *)calloc((size_t)d2->ncc, sizeof(int8_t));
  d2->dots[0] = 1;
  assert(CannedCobordismImpl_check(d2));

  /* Compose: d2 on top of d1 (d1.top == capA == d2.bottom) */
  CannedCobordism *composed = CannedCobordismImpl_compose_vertical(d1, d2);
  assert(composed != NULL);

  CannedCobordismImplData *rd = (CannedCobordismImplData *)composed->impl_data;
  assert(CannedCobordismImpl_check(rd) == true);
  /* Result: top = d2.top = capB,  bottom = d1.bottom = capB */
  assert(Cap_equals(rd->top, capB));
  assert(Cap_equals(rd->bottom, capB));

  CannedCobordism_free(composed);
  CannedCobordismImpl_free(d1);
  CannedCobordismImpl_free(d2);
  Cap_free(capA);
  Cap_free(capB);
  printf("  [PASS] test_compose_vertical_nontrivial\n");
}

/* ================================================================
 *  Test 18: compose_horizontal
 *
 *  Take two 4-edge isomorphism cobordisms, compose horizontally
 *  joining 2 edges. Result should have n = 4 + 4 - 2*2 = 4 edges.
 * ================================================================ */
static void test_compose_horizontal(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);

  CannedCobordism *iso1 = CannedCobordismImpl_isomorphism(c);
  CannedCobordism *iso2 = CannedCobordismImpl_isomorphism(c);

  CannedCobordismImplData *d1 = (CannedCobordismImplData *)iso1->impl_data;
  CannedCobordismImplData *d2 = (CannedCobordismImplData *)iso2->impl_data;

  CannedCobordism *composed =
      CannedCobordismImpl_compose_horizontal(d1, 0, d2, 0, 2);
  assert(composed != NULL);

  CannedCobordismImplData *rd = (CannedCobordismImplData *)composed->impl_data;
  assert(rd->n == 4); /* 4 + 4 - 2*2 */
  assert(CannedCobordismImpl_check(rd) == true);

  CannedCobordism_free(composed);
  CannedCobordism_free(iso1);
  CannedCobordism_free(iso2);
  Cap_free(c);
  printf("  [PASS] test_compose_horizontal\n");
}

/* ================================================================
 *  Test 19: as_CannedCobordism — vtable dispatch
 *
 *  Wrap an impl, then call through the polymorphic interface.
 * ================================================================ */
static void test_vtable_dispatch(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);

  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  assert(iso != NULL);

  /* isIsomorphism through vtable */
  assert(CannedCobordism_isIsomorphism(iso) == true);

  /* source / target through vtable */
  assert(Cap_equals(CannedCobordism_source(iso), c));
  assert(Cap_equals(CannedCobordism_target(iso), c));

  /* reverseMaps through vtable */
  CannedCobordism_reverseMaps(iso);
  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;
  assert(d->reverse_maps_done == true);

  /* compose through vtable (identity ∘ identity) */
  CannedCobordism *composed = CannedCobordism_compose(iso, iso);
  assert(composed != NULL);
  assert(CannedCobordism_isIsomorphism(composed) == true);

  CannedCobordism_free(composed);
  CannedCobordism_free(iso);
  Cap_free(c);
  printf("  [PASS] test_vtable_dispatch\n");
}

/* ================================================================
 *  Test 20: free — no crash on populated & NULL
 * ================================================================ */
static void test_free_no_crash(void) {
  int p[] = {1, 0, 3, 2};
  Cap *c = make_cap(4, 0, p);

  CannedCobordism *iso = CannedCobordismImpl_isomorphism(c);
  CannedCobordismImplData *d = (CannedCobordismImplData *)iso->impl_data;

  /* Populate reverse maps before freeing */
  CannedCobordismImpl_reverseMaps(d);
  CannedCobordism_free(iso);

  /* Free NULL should be safe */
  CannedCobordismImpl_free(NULL);
  CannedCobordism_free(NULL);

  Cap_free(c);
  printf("  [PASS] test_free_no_crash\n");
}

/* ================================================================
 *  main
 * ================================================================ */
int main(void) {
  printf("Running CannedCobordismImpl tests...\n\n");

  test_create_simple();
  test_create_different_caps();
  test_isomorphism_factory();
  test_isIsomorphism_false();
  test_check_valid();
  test_check_invalid();
  test_reverseMaps();
  test_reverseMaps_idempotent();
  test_equals_self();
  test_equals_cloned();
  test_equals_different();
  test_hashCode_same();
  test_hashCode_caching();
  test_compareTo_self();
  test_compareTo_ordering();
  test_compose_vertical_identity();
  test_compose_vertical_nontrivial();
  test_compose_horizontal();
  test_vtable_dispatch();
  test_free_no_crash();

  printf("\nAll CannedCobordismImpl tests passed!\n");
  return 0;
}
