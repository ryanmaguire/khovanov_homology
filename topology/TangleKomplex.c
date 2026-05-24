#include "TangleKomplex.h"
#include "CannedCobordismImpl.h"
#include "LCCC.h"
#include <stdlib.h>
/* Helper to horizontally compose two Caps */
static Cap *Cap_compose_tangle(const Cap *a, const Cap *b, int n_strands) {
  /* Join the n_strands right endpoints of a with the n_strands left endpoints of b.
     In our convention, endpoints 0..n-1 are left, n..2n-1 are right. */
  return Cap_compose(a, n_strands, b, 0, n_strands, NULL);
}
/* Helper to horizontally compose two SmoothingColumns */
static SmoothingColumn *SmoothingColumn_tensorProduct(SmoothingColumn *A, SmoothingColumn *B, int n_strands) {
  SmoothingColumn *C = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  C->n = A->n * B->n;
  if (C->n > 0) {
    C->numbers = (int *)malloc((size_t)C->n * sizeof(int));
    C->smoothings = (Cap **)malloc((size_t)C->n * sizeof(Cap *));
  } else {
    C->numbers = NULL;
    C->smoothings = NULL;
  }
  int k = 0;
  for (int i = 0; i < A->n; i++) {
    for (int j = 0; j < B->n; j++) {
      C->numbers[k] = A->numbers[i] + B->numbers[j];
      /* The tensor product Cap */
      Cap *z = Cap_compose_tangle(A->smoothings[i], B->smoothings[j], n_strands);
      
      /* Bar-Natan algorithm: if new cycles are formed during horizontal composition,
         they must be delooped. However, the basic Cap_compose just increments ncycles.
         The delooping will happen during Komplex_reduce_with_oracle because the cycles
         become isolated circles. Actually, in Bar-Natan's formal complex, circles are 
         delooped into two basis elements (+1 and -1 shift).

      C->smoothings[k] = z;
      k++;
    }
  }
  return C;
}
/* Helper to compute d_L \otimes id_R */
static CobMatrix *CobMatrix_tensorProductMap1(CobMatrix *d_left, SmoothingColumn *col_right, int n_strands) {
  SmoothingColumn *new_src = SmoothingColumn_tensorProduct(d_left->source, col_right, n_strands);
  SmoothingColumn *new_tgt = SmoothingColumn_tensorProduct(d_left->target, col_right, n_strands);
  
  CobMatrix *res = CobMatrix_create(new_src, new_tgt, true);
  
  for (int i = 0; i < d_left->target->n; i++) {
    MatrixEntry *entry = d_left->entries[i].head;
    while (entry != NULL) {
      int c = entry->column_index;
      
      /* For each right basis element j, we add the entry f \otimes id */
      for (int j = 0; j < col_right->n; j++) {
        int row_idx = i * col_right->n + j;
        int col_idx = c * col_right->n + j;
        
        /* Create f \otimes id */
        LCCC *f_tensor_id = LCCC_createZero();
        LCCCTerm *term = entry->value->head;
        while (term != NULL) {
          CannedCobordism *f = term->cobordism;
          CannedCobordism *id = CannedCobordismImpl_isomorphism(col_right->smoothings[j]);
          
          CannedCobordism *ftid = f->compose_partial(f, n_strands, id, 0, n_strands);
          
          /* Add to LCCC */
          LCCC *single = LCCC_createSingle(ftid);
          LCCC *sum = LCCC_add(f_tensor_id, single);
          LCCC_free(f_tensor_id);
          LCCC_free(single);
          f_tensor_id = sum;
          
          id->free_impl(id);
          term = term->next;
        }
        
        CobMatrix_putEntry(res, row_idx, col_idx, f_tensor_id);
      }
      entry = entry->next;
    }
  }
  return res;
}
/* Helper to compute (-1)^{deg_left} id_L \otimes d_R */
static CobMatrix *CobMatrix_tensorProductMap2(SmoothingColumn *col_left, CobMatrix *d_right, int n_strands, int deg_left) {
  SmoothingColumn *new_src = SmoothingColumn_tensorProduct(col_left, d_right->source, n_strands);
  SmoothingColumn *new_tgt = SmoothingColumn_tensorProduct(col_left, d_right->target, n_strands);
  
  CobMatrix *res = CobMatrix_create(new_src, new_tgt, true);
  
  /* Over F_2, the sign (-1)^deg_left is irrelevant. */
  (void)deg_left;
  
  for (int i = 0; i < d_right->target->n; i++) {
    MatrixEntry *entry = d_right->entries[i].head;
    while (entry != NULL) {
      int c = entry->column_index;
      
      /* For each left basis element j, we add id \otimes g */
      for (int j = 0; j < col_left->n; j++) {
        int row_idx = j * d_right->target->n + i;
        int col_idx = j * d_right->source->n + c;
        
        LCCC *id_tensor_g = LCCC_createZero();
        LCCCTerm *term = entry->value->head;
        while (term != NULL) {
          CannedCobordism *g = term->cobordism;
          CannedCobordism *id = CannedCobordismImpl_isomorphism(col_left->smoothings[j]);
          
          CannedCobordism *idtg = id->compose_partial(id, n_strands, g, 0, n_strands);
          
          LCCC *single = LCCC_createSingle(idtg);
          LCCC *sum = LCCC_add(id_tensor_g, single);
          LCCC_free(id_tensor_g);
          LCCC_free(single);
          id_tensor_g = sum;
          
          id->free_impl(id);
          term = term->next;
        }
        
        CobMatrix_putEntry(res, row_idx, col_idx, id_tensor_g);
      }
      entry = entry->next;
    }
  }
  return res;
}
Komplex *Komplex_compose_tangles(Komplex *L, Komplex *R, int n_strands) {
  int new_length = L->length + R->length - 1;
  Komplex *T = Komplex_create(new_length);
  
  /* Create chain groups T_k = \bigoplus_{i+j=k} (L_i \otimes R_j) */
  for (int k = 0; k < new_length; k++) {
    /* Count total size */
    int total_n = 0;
    for (int i = 0; i < L->length; i++) {
      int j = k - i;
      if (j >= 0 && j < R->length) {
        total_n += L->chain_groups[i]->n * R->chain_groups[j]->n;
      }
    }
    
    SmoothingColumn *Ck = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
    Ck->n = total_n;
    if (total_n > 0) {
      Ck->numbers = (int *)malloc((size_t)total_n * sizeof(int));
      Ck->smoothings = (Cap **)malloc((size_t)total_n * sizeof(Cap *));
    } else {
      Ck->numbers = NULL;
      Ck->smoothings = NULL;
    }
    
    int idx = 0;
    for (int i = 0; i < L->length; i++) {
      int j = k - i;
      if (j >= 0 && j < R->length) {
        SmoothingColumn *Li = L->chain_groups[i];
        SmoothingColumn *Rj = R->chain_groups[j];
        for (int a = 0; a < Li->n; a++) {
          for (int b = 0; b < Rj->n; b++) {
            Ck->numbers[idx] = Li->numbers[a] + Rj->numbers[b];
            Ck->smoothings[idx] = Cap_compose_tangle(Li->smoothings[a], Rj->smoothings[b], n_strands);
            idx++;
          }
        }
      }
    }
    T->chain_groups[k] = Ck;
  }
  
  /* Create differentials D_k : T_k -> T_{k+1} */
  for (int k = 0; k < new_length - 1; k++) {
    CobMatrix *Dk = CobMatrix_create(T->chain_groups[k], T->chain_groups[k+1], true);
    
    int row_offset = 0;
    for (int i_out = 0; i_out < L->length; i_out++) {
      int j_out = (k + 1) - i_out;
      if (j_out >= 0 && j_out < R->length) {
        
        int col_offset = 0;
        for (int i_in = 0; i_in < L->length; i_in++) {
          int j_in = k - i_in;
          if (j_in >= 0 && j_in < R->length) {
            
            /* D_k restricts to d_L \otimes id + id \otimes d_R */
            /* d_L part: if i_out == i_in + 1 and j_out == j_in */
            if (i_out == i_in + 1 && j_out == j_in) {
              CobMatrix *dL = L->differentials[i_in];
              CobMatrix *d_tensor_id = CobMatrix_tensorProductMap1(dL, R->chain_groups[j_in], n_strands);
              
              /* Add block to Dk */
              for (int r = 0; r < d_tensor_id->target->n; r++) {
                MatrixEntry *entry = d_tensor_id->entries[r].head;
                while (entry != NULL) {
                  CobMatrix_addEntry(Dk, row_offset + r, col_offset + entry->column_index, LCCC_clone(entry->value));
                  entry = entry->next;
                }
              }
              CobMatrix_free(d_tensor_id);
            }
            
            /* d_R part: if i_out == i_in and j_out == j_in + 1 */
            if (i_out == i_in && j_out == j_in + 1) {
              CobMatrix *dR = R->differentials[j_in];
              CobMatrix *id_tensor_d = CobMatrix_tensorProductMap2(L->chain_groups[i_in], dR, n_strands, i_in);
              
              /* Add block to Dk */
              for (int r = 0; r < id_tensor_d->target->n; r++) {
                MatrixEntry *entry = id_tensor_d->entries[r].head;
                while (entry != NULL) {
                  CobMatrix_addEntry(Dk, row_offset + r, col_offset + entry->column_index, LCCC_clone(entry->value));
                  entry = entry->next;
                }
              }
              CobMatrix_free(id_tensor_d);
            }
            
            col_offset += L->chain_groups[i_in]->n * R->chain_groups[j_in]->n;
          }
        }
        
        row_offset += L->chain_groups[i_out]->n * R->chain_groups[j_out]->n;
      }
    }
    
    T->differentials[k] = Dk;
  }
  
  return T;
}
Komplex *Komplex_identityBraid(int n_strands) {
  Komplex *k = Komplex_create(1);
  SmoothingColumn *c0 = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  c0->n = 1;
  c0->numbers = (int *)malloc(sizeof(int));
  c0->numbers[0] = 0;
  
  c0->smoothings = (Cap **)malloc(sizeof(Cap *));
  Cap *id_cap = Cap_create(2 * n_strands, 0);
  for (int i = 0; i < n_strands; i++) {
    id_cap->pairings[i] = i + n_strands;
    id_cap->pairings[i + n_strands] = i;
  }
  c0->smoothings[0] = id_cap;
  
  k->chain_groups[0] = c0;
  return k;
}
Komplex *Komplex_singleCrossing(int n_strands, int crossing_index, bool positive) {
  Komplex *k = Komplex_create(2);
  
  /* The crossing occurs between strands: crossing_index and crossing_index + 1 */
  /* We build the 0-smoothing and 1-smoothing on 2*n_strands points */
  
  Cap *cap0 = Cap_create(2 * n_strands, 0);
  Cap *cap1 = Cap_create(2 * n_strands, 0);
  
  for (int i = 0; i < n_strands; i++) {
    if (i == crossing_index || i == crossing_index + 1) {
      /* 0-smoothing: identity matching */
      cap0->pairings[i] = i + n_strands;
      cap0->pairings[i + n_strands] = i;
      
      /* 1-smoothing: horizontal matching */
      int other = (i == crossing_index) ? crossing_index + 1 : crossing_index;
      cap1->pairings[i] = other;
      cap1->pairings[i + n_strands] = other + n_strands;
    } else {
      /* Identity for untouched strands */
      cap0->pairings[i] = i + n_strands;
      cap0->pairings[i + n_strands] = i;
      cap1->pairings[i] = i + n_strands;
      cap1->pairings[i + n_strands] = i;
    }
  }
  
  SmoothingColumn *c0 = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  c0->n = 1;
  c0->numbers = (int *)malloc(sizeof(int));
  c0->smoothings = (Cap **)malloc(sizeof(Cap *));
  c0->smoothings[0] = positive ? cap0 : cap1;
  c0->numbers[0] = positive ? 0 : -1; /* Example shift, depends on FPT convention */
  
  SmoothingColumn *c1 = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  c1->n = 1;
  c1->numbers = (int *)malloc(sizeof(int));
  c1->smoothings = (Cap **)malloc(sizeof(Cap *));
  c1->smoothings[0] = positive ? cap1 : cap0;
  c1->numbers[0] = positive ? 1 : 0;
  
  k->chain_groups[0] = c0;
  k->chain_groups[1] = c1;
  
  /* Saddle cobordism */
  CannedCobordismImplData *impl = CannedCobordismImpl_create(c0->smoothings[0], c1->smoothings[0]);
  /* Initialize simple CCs. Isomorphism and saddle don't have dots/genus initially */
  impl->ncc = impl->nbc; 
     
  impl->ncc = impl->nbc;
  for (int i = 0; i < impl->nbc; i++) impl->connectedComponent[i] = i;
  impl->dots = (int8_t *)calloc((size_t)impl->ncc, sizeof(int8_t));
  impl->genus = (int8_t *)calloc((size_t)impl->ncc, sizeof(int8_t));
  
  CannedCobordism *saddle = CannedCobordismImpl_as_CannedCobordism(impl);
  
  CobMatrix *d = CobMatrix_create(c0, c1, true);
  CobMatrix_putEntry(d, 0, 0, LCCC_createSingle(saddle));
  k->differentials[0] = d;
  
  return k;
}
