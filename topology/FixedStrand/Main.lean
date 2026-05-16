/-
  Fixed-Strand Conjecture for 4-strand braids: Main theorem
  
  For any 4-strand braid b, the independence complex of the augmented rhomboid graph
  collapses to a connected complex with no free pairs.
-/
import FixedStrand.Defs
import FixedStrand.Collapse
import FixedStrand.Connectivity
import FixedStrand.SmallCases
import FixedStrand.LargeCases

set_option maxHeartbeats 800000

/-! ## Main connectivity theorem -/

/-- The independence complex of the augmented rhomboid graph is connected for any 4-strand braid. -/
theorem braid_4strand_connected (b : Braid 4) :
    is_connected (IndependenceComplex (AugmentedRhomboidGraph b)) := by
  by_cases h : b.length < 3
  · exact small_braids_connected b h
  · exact large_braids_connected b (by omega)

/-! ## Faces of independence complexes are finite -/

lemma indep_faces_finite (b : Braid 4) :
    (IndependenceComplex (AugmentedRhomboidGraph b)).faces.Finite := by
  apply Set.Finite.subset (Set.finite_univ.subset (Set.subset_univ _))
  intro s _; exact Set.mem_univ s

/-! ## Main theorem -/

/-- Fixed-Strand Conjecture for 4-strand braids:
    For any 4-strand braid b, the independence complex of its augmented rhomboid graph
    collapses to a connected complex with no free pairs.
    
    This implies the independence complex has the homotopy type of a wedge of spheres,
    which is the key claim of the Fixed-Strand Conjecture. -/
theorem fixed_strand_conjecture_4 (b : Braid 4) :
    ∃ L, Collapses (IndependenceComplex (AugmentedRhomboidGraph b)) L ∧
         is_connected L ∧ has_no_free_pairs L := by
  let K := IndependenceComplex (AugmentedRhomboidGraph b)
  refine ⟨collapse_all K, collapse_all_collapses K, ?_, ?_⟩
  · exact Collapses_preserves_connected (collapse_all_collapses K) (braid_4strand_connected b)
  · exact collapse_all_has_no_free_pairs K (indep_faces_finite b)
