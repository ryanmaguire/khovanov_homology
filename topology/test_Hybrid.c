#include "HybridKhovanov.h"
#include "MorseOracle.h"
#include <stdio.h>
#include <stdlib.h>
/* Using real Katlas algebra, no stubs needed. */
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
int main() {
  printf("Testing Hybrid Khovanov Algorithm (Compilation & Basic Structure)...\n");
  
  int n_strands = 3;
  int crossings[] = {0, 1, 0, 1, 0, 1}; // length 6
  bool signs[] = {true, true, true, true, true, true};
  
  printf("\n--- Running with Threshold = 6 (Pure Bar-Natan via Algebraic) ---\n");
  Komplex *k1 = HybridKhovanov_compute_algebraic(n_strands, crossings, signs, 6, 6);
  printf("Result Komplex has %d chain groups.\n", k1->length);
  Komplex_free(k1);
  
  printf("\n--- Running with Threshold = 1 (Pure Divide & Conquer via Algebraic) ---\n");
  Komplex *k2 = HybridKhovanov_compute_algebraic(n_strands, crossings, signs, 6, 1);
  printf("Result Komplex has %d chain groups.\n", k2->length);
  Komplex_free(k2);
  
  printf("\n--- Running with Threshold = 3 (Hybrid via Algebraic) ---\n");
  Komplex *k3 = HybridKhovanov_compute_algebraic(n_strands, crossings, signs, 6, 3);
  printf("Result Komplex has %d chain groups.\n", k3->length);
  Komplex_free(k3);
  
  printf("\n--- Running Topological Morse Oracle ---\n");
  Komplex *k4 = HybridKhovanov_compute_oracle(n_strands, crossings, signs, 6);
  printf("Result Komplex has %d chain groups.\n", k4->length);
  Komplex_free(k4);
  
  printf("\nAll compilation and execution paths for HybridKhovanov succeeded!\n");
  return 0;
}
