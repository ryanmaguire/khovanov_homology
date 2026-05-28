#include "HybridKhovanov.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
/* Stub Cap / SmoothingColumn deep copy for the matrix stubs */
SmoothingColumn *SmoothingColumn_clone(SmoothingColumn *col) {
  SmoothingColumn *clone = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  clone->n = col->n;
  clone->numbers = (int *)malloc((size_t)col->n * sizeof(int));
  clone->smoothings = (Cap **)malloc((size_t)col->n * sizeof(Cap *));
  return clone;
}
bool SmoothingColumn_equals(SmoothingColumn *a, SmoothingColumn *b) {
  return a->n == b->n;
}
int main(int argc, char **argv) {
  int length = 8;
  if (argc > 1) {
    length = atoi(argv[1]);
  }
  int n_strands = 3;
  int *crossings = (int*)malloc(length * sizeof(int));
  bool *signs = (bool*)malloc(length * sizeof(bool));
  
  for (int i = 0; i < length; i++) {
    crossings[i] = i % (n_strands - 1);
    signs[i] = (i % 2 == 0);
  }
  printf("=== Benchmarking Fast Khovanov Algorithms ===\n");
  printf("Braid: %d strands, %d crossings\n", n_strands, length);
  clock_t start, end;
  double time_used;
  /* 1. Pure Bar-Natan (Threshold = Length) */
  printf("\n1. Pure Bar-Natan (Threshold = %d)\n", length);
  start = clock();
  Komplex *k_bn = HybridKhovanov_compute_algebraic(n_strands, crossings, signs, length, length);
  end = clock();
  time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("   Time: %f seconds\n", time_used);
  Komplex_free(k_bn);
  /* 2. Pure Divide & Conquer (Threshold = 1) */
  printf("\n2. Pure Divide & Conquer (Threshold = 1)\n");
  start = clock();
  Komplex *k_dc = HybridKhovanov_compute_algebraic(n_strands, crossings, signs, length, 1);
  end = clock();
  time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("   Time: %f seconds\n", time_used);
  Komplex_free(k_dc);
  /* 3. Hybrid (Threshold = 4) */
  int thresh = (length > 4) ? 4 : 2;
  printf("\n3. Hybrid D&C (Threshold = %d)\n", thresh);
  start = clock();
  Komplex *k_hy = HybridKhovanov_compute_algebraic(n_strands, crossings, signs, length, thresh);
  end = clock();
  time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("   Time: %f seconds\n", time_used);
  Komplex_free(k_hy);
  free(crossings);
  free(signs);
  return 0;
}
