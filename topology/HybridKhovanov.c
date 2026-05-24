#include "HybridKhovanov.h"
#include "TangleKomplex.h"
#include "MorseOracle.h"
#include <stdlib.h>
Komplex *HybridKhovanov_compute_algebraic(int n_strands, int *crossings, bool *signs, int length, int threshold) {
  if (length == 0) {
    return Komplex_identityBraid(n_strands);
  }
  
  if (length <= threshold) {
    /* Base Case: Bar-Natan iterative scan with greedy reduction */
    Komplex *K = Komplex_identityBraid(n_strands);
    
    for (int i = 0; i < length; i++) {
      Komplex *K_cross = Komplex_singleCrossing(n_strands, crossings[i], signs[i]);
      Komplex *K_next = Komplex_compose_tangles(K, K_cross, n_strands);
      
      Komplex_free(K);
      Komplex_free(K_cross);
      
      K = K_next;
      
      /* Apply Bar-Natan's greedy algebraic reduction step */
      Komplex_greedyReduce(K);
    }
    
    return K;
  } else {
    /* Recursive Step: Divide and Conquer with greedy reduction */
    int mid = length / 2;
    
    Komplex *K_L = HybridKhovanov_compute_algebraic(n_strands, crossings, signs, mid, threshold);
    Komplex *K_R = HybridKhovanov_compute_algebraic(n_strands, crossings + mid, signs + mid, length - mid, threshold);
    
    Komplex *K_total = Komplex_compose_tangles(K_L, K_R, n_strands);
    
    Komplex_free(K_L);
    Komplex_free(K_R);
    
    /* Apply greedy algebraic reduction on the combined tangle */
    Komplex_greedyReduce(K_total);
       
    return K_total;
  }
}
Komplex *HybridKhovanov_compute_oracle(int n_strands, int *crossings, bool *signs, int length) {
  /* 1. Build the global unreduced chain complex */
  Komplex *K = Komplex_identityBraid(n_strands);
  for (int i = 0; i < length; i++) {
    Komplex *K_cross = Komplex_singleCrossing(n_strands, crossings[i], signs[i]);
    Komplex *K_next = Komplex_compose_tangles(K, K_cross, n_strands);
    Komplex_free(K);
    Komplex_free(K_cross);
    K = K_next;
  }
  
  if (length > 0) {
    /* 2. Build the Braid object for the topological Oracle */
    /* Note: We currently assume signs are positive for the topological graph 
       as per MorseOracle.h limitations, but pass the generators. */
    Braid *b = Braid_create(n_strands, length, crossings);
    
    /* 3. Compute schedule and reduce */
    CollapseSchedule *schedule = MorseOracle_fullSchedule(b);
    Komplex_reduce_with_oracle(K, schedule);
    
    CollapseSchedule_free(schedule);
    Braid_free(b);
  }
  
  return K;
}
