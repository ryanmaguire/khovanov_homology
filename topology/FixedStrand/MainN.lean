/-
  Fixed-Strand Conjecture: Combined theorem

  This file combines:
  1. The complete 4-strand proof (from FixedStrand.Main)
  2. The general degree bound (from FixedStrandGeneral.DegreeBound)
  3. The sparsity-based connectivity argument (from FixedStrand.LargeCases)

  to produce:
  - `fixed_strand_conjecture_general`: For any n-strand braid with n*(L+1) > 2*6+2,
    the independence complex collapses to a connected complex with no free pairs.
  - Corollaries for specific n values.
  - A unified `fixed_strand_conjecture_4` that combines small cases with the general argument.
-/
import FixedStrandGeneral.DegreeBound
import FixedStrand.Collapse
import FixedStrand.Connectivity
import FixedStrand.LargeCases
import FixedStrand.SmallCases
import FixedStrand.Main

set_option maxHeartbeats 800000

open scoped Classical

/-! ## Vertex count -/

/-- The number of vertices in the augmented rhomboid graph is n * (b.length + 1). -/
lemma aug_vertex_count (n : ℕ) (b : Braid n) :
    Fintype.card (Fin n × Fin (b.length + 1)) = n * (b.length + 1) := by
  simp [Fintype.card_prod, Fintype.card_fin]

/-! ## General connectivity theorem -/

/-- For any n-strand braid with n * (b.length + 1) > 2 * 6 + 2 = 14, the independence
    complex of the augmented rhomboid graph is connected.

    This uses: degree ≤ 6 (general n) + compl_connected_of_low_degree. -/
theorem general_braids_connected (n : ℕ) (b : Braid n)
    (hsize : n * (b.length + 1) > 2 * 6 + 2) :
    is_connected (IndependenceComplex (AugmentedRhomboidGraph b)) := by
  apply indep_connected_of_compl_connected
  apply compl_connected_of_low_degree (k := 6)
  · intro v; exact aug_degree_le_6_general n b v
  · rw [aug_vertex_count]; exact hsize

/-! ## Faces finiteness for general n -/

lemma indep_faces_finite_general (n : ℕ) (b : Braid n) :
    (IndependenceComplex (AugmentedRhomboidGraph b)).faces.Finite := by
  apply Set.Finite.subset (Set.finite_univ.subset (Set.subset_univ _))
  intro s _; exact Set.mem_univ s

/-! ## Main theorem: Fixed-Strand Conjecture for general n -/

/-- Fixed-Strand Conjecture for general n-strand braids:
    For any n-strand braid b with n * (b.length + 1) > 14, the independence complex
    of its augmented rhomboid graph collapses to a connected complex with no free pairs.

    The condition n * (b.length + 1) > 14 is satisfied whenever the braid has enough
    vertices for the sparsity argument to apply. In particular:
    - For n = 4: holds for b.length ≥ 3 (since 4 * 4 = 16 > 14)
    - For n = 5: holds for b.length ≥ 2 (since 5 * 3 = 15 > 14)
    - For n = 8: holds for b.length ≥ 1 (since 8 * 2 = 16 > 14)
    - For n ≥ 15: holds for ALL braids   (since 15 * 1 = 15 > 14) -/
theorem fixed_strand_conjecture_general (n : ℕ) (b : Braid n)
    (hsize : n * (b.length + 1) > 14) :
    ∃ L, Collapses (IndependenceComplex (AugmentedRhomboidGraph b)) L ∧
         is_connected L ∧ has_no_free_pairs L := by
  let K := IndependenceComplex (AugmentedRhomboidGraph b)
  refine ⟨collapse_all K, collapse_all_collapses K, ?_, ?_⟩
  · exact Collapses_preserves_connected (collapse_all_collapses K)
      (general_braids_connected n b hsize)
  · exact collapse_all_has_no_free_pairs K (indep_faces_finite_general n b)

/-! ## Corollaries for specific n values -/

/-- For n ≥ 15, the Fixed-Strand Conjecture holds for ALL braids (including empty). -/
theorem fixed_strand_conjecture_large_n (n : ℕ) (hn : n ≥ 15) (b : Braid n) :
    ∃ L, Collapses (IndependenceComplex (AugmentedRhomboidGraph b)) L ∧
         is_connected L ∧ has_no_free_pairs L := by
  apply fixed_strand_conjecture_general
  have : b.length + 1 ≥ 1 := by omega
  nlinarith

/-- For n ≥ 5 and braids of length ≥ 2, the conjecture holds. -/
theorem fixed_strand_conjecture_n_ge_5 (n : ℕ) (hn : n ≥ 5) (b : Braid n)
    (hlen : b.length ≥ 2) :
    ∃ L, Collapses (IndependenceComplex (AugmentedRhomboidGraph b)) L ∧
         is_connected L ∧ has_no_free_pairs L := by
  apply fixed_strand_conjecture_general
  nlinarith

/-- For n ≥ 4 and braids of length ≥ 3, the conjecture holds.
    This recovers the long-braid part of the 4-strand result via pure sparsity. -/
theorem fixed_strand_conjecture_n_ge_4_long (n : ℕ) (hn : n ≥ 4) (b : Braid n)
    (hlen : b.length ≥ 3) :
    ∃ L, Collapses (IndependenceComplex (AugmentedRhomboidGraph b)) L ∧
         is_connected L ∧ has_no_free_pairs L := by
  apply fixed_strand_conjecture_general
  nlinarith

/-! ## The complete 4-strand theorem (unified) -/

/-- The Fixed-Strand Conjecture for 4-strand braids, unified from all cases.

    This is exactly `fixed_strand_conjecture_4` from FixedStrand.Main, re-exported
    here so that all results are accessible from one module. The 4-strand proof
    combines:
    - Small cases (length ≤ 2): BFS + native_decide via SmallCases.lean
    - Large cases (length ≥ 3): degree bound + sparsity via LargeCases.lean

    The general degree bound in this module shows the same approach works for
    arbitrary n, with only the small-case threshold changing. -/
theorem fixed_strand_conjecture_4_complete (b : Braid 4) :
    ∃ L, Collapses (IndependenceComplex (AugmentedRhomboidGraph b)) L ∧
         is_connected L ∧ has_no_free_pairs L :=
  fixed_strand_conjecture_4 b
