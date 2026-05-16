/-
  Fixed-Strand Conjecture: General degree bound for n-strand braids

  Key result: Every vertex in the augmented rhomboid graph of an n-strand braid
  has degree at most 6, regardless of n.

  Proof strategy: Neighbors of v = (s, t) come from three disjoint time levels:
  - Time t (vertical edges): at most 2 neighbors (strands s±1)
  - Time t+1 (horizontal + crossing forward): at most 2 neighbors
  - Time t-1 (horizontal + crossing backward): at most 2 neighbors
  Total: at most 6.
-/
import FixedStrand.Defs

set_option maxHeartbeats 800000

open scoped Classical

/-! ## Adjacency implies proximity -/

/-- Adjacency in the augmented rhomboid graph requires time coordinates to differ
    by at most 1. -/
lemma aug_adj_time_close (n : ℕ) (b : Braid n)
    (v w : Fin n × Fin (b.length + 1))
    (hadj : (AugmentedRhomboidGraph b).Adj v w) :
    v.2 = w.2 ∨ v.2.val + 1 = w.2.val ∨ w.2.val + 1 = v.2.val := by
  cases hadj
  unfold AugRhomboidAdj at *; aesop

/-- Adjacency requires strand coordinates to differ by at most 1. -/
lemma aug_adj_strand_close (n : ℕ) (b : Braid n)
    (v w : Fin n × Fin (b.length + 1))
    (hadj : (AugmentedRhomboidGraph b).Adj v w) :
    v.1 = w.1 ∨ v.1.val + 1 = w.1.val ∨ w.1.val + 1 = v.1.val := by
  unfold AugmentedRhomboidGraph at hadj
  unfold AugRhomboidAdj at hadj
  grind +suggestions

/-! ## At most 2 neighbors at each time level -/

/-- At most 2 neighbors at the same time step (from vertical edges only). -/
lemma card_neighbors_same_time (n : ℕ) (b : Braid n)
    (v : Fin n × Fin (b.length + 1)) :
    (Finset.univ.filter fun w => (AugmentedRhomboidGraph b).Adj v w ∧ w.2 = v.2).card ≤ 2 := by
  refine' le_trans ( Finset.card_le_card _ ) _
  exact { ( ⟨ if v.1.val = 0 then 0 else v.1.val - 1, by
    split_ifs <;> omega ⟩, v.2 ), ( ⟨ if v.1.val = n - 1 then n - 1 else v.1.val + 1, by
    grind ⟩, v.2 ) }
  · intro w hw; simp_all +decide [ AugmentedRhomboidGraph ]
    rcases hw with ⟨ ⟨ hne, h | h ⟩, hw ⟩ <;> simp_all +decide [ AugRhomboidAdj ]
    · grind
    · grind
  · exact Finset.card_insert_le _ _

/-- At most 2 neighbors at the next time step (from horizontal + crossing forward). -/
lemma card_neighbors_next_time (n : ℕ) (b : Braid n)
    (v : Fin n × Fin (b.length + 1)) :
    (Finset.univ.filter fun w => (AugmentedRhomboidGraph b).Adj v w ∧ v.2.val + 1 = w.2.val).card ≤ 2 := by
  rcases v with ⟨ i, j ⟩
  by_cases h : j.val < b.length
  · refine' le_trans ( Finset.card_le_card _ ) _
    exact { ( i, ⟨ j.val + 1, by linarith ⟩ ) } ∪ ( if h : ( b.get ⟨ j.val, h ⟩ ).val = i.val then { ( ⟨ i.val + 1, by
      grind ⟩, ⟨ j.val + 1, by linarith ⟩ ) } else ∅ ) ∪ ( if h : ( b.get ⟨ j.val, h ⟩ ).val = i.val - 1 ∧ i.val > 0 then { ( ⟨ i.val - 1, by
      exact lt_of_le_of_lt ( Nat.pred_le _ ) i.2 ⟩, ⟨ j.val + 1, by linarith ⟩ ) } else ∅ )
    · all_goals generalize_proofs at *
      intro w hw; simp_all +decide [ AugmentedRhomboidGraph ]
      rcases hw with ⟨ ⟨ hw₁, hw₂ ⟩, hw₃ ⟩ ; rcases hw₂ with ( hw₂ | hw₂ ) <;> simp_all +decide [ AugRhomboidAdj ]
      · grind
      · grind
    · grind +splitIndPred
  · simp_all +decide [ Fin.eq_last_of_not_lt h ]
    rw [ Finset.card_eq_zero.mpr ] <;> norm_num
    exact fun a b h => ne_of_gt ( Nat.lt_succ_of_le ( Fin.is_le b ) )

/-- At most 2 neighbors at the previous time step (from horizontal + crossing backward). -/
lemma card_neighbors_prev_time (n : ℕ) (b : Braid n)
    (v : Fin n × Fin (b.length + 1)) :
    (Finset.univ.filter fun w => (AugmentedRhomboidGraph b).Adj v w ∧ w.2.val + 1 = v.2.val).card ≤ 2 := by
  unfold AugmentedRhomboidGraph
  rcases v with ⟨ x, y ⟩
  rcases y with ⟨ _ | y, hy ⟩
  · exact Finset.card_eq_zero.mpr ( by aesop ) |> fun h => h.le.trans ( by norm_num )
  · refine' le_trans ( Finset.card_le_card _ ) _
    exact { ( x, ⟨ y, by linarith ⟩ ), ( if x.val = ( b.get ⟨ y, by linarith ⟩ ).val + 1 then ⟨ ( b.get ⟨ y, by linarith ⟩ ).val, by
      grind ⟩ else if x.val = ( b.get ⟨ y, by linarith ⟩ ).val then ⟨ ( b.get ⟨ y, by linarith ⟩ ).val + 1, by
      grind ⟩ else x, ⟨ y, by linarith ⟩ ) }
    · simp +decide [ Finset.subset_iff, AugRhomboidAdj ]
      rintro a b h₁ h₂ h₃; rcases b with ⟨ _ | b, hb ⟩ <;> simp_all +decide [ Fin.ext_iff ]
      · grind
      · grind
    · exact Finset.card_insert_le _ _

/-! ## Main degree bound -/

/-- The degree of any vertex in the augmented rhomboid graph is at most 6,
    regardless of the number of strands n. -/
theorem aug_degree_le_6_general (n : ℕ) (b : Braid n)
    (v : Fin n × Fin (b.length + 1)) :
    (Finset.univ.filter fun w => (AugmentedRhomboidGraph b).Adj v w).card ≤ 6 := by
  -- Partition the neighbor set by time coordinate into three disjoint subsets
  have h_partition : Finset.univ.filter (fun w => (AugmentedRhomboidGraph b).Adj v w) ⊆
    (Finset.univ.filter (fun w => (AugmentedRhomboidGraph b).Adj v w ∧ w.2 = v.2)) ∪
      (Finset.univ.filter (fun w => (AugmentedRhomboidGraph b).Adj v w ∧ v.2.val + 1 = w.2.val)) ∪
      (Finset.univ.filter (fun w => (AugmentedRhomboidGraph b).Adj v w ∧ w.2.val + 1 = v.2.val)) := by
    intro w hw; have := aug_adj_time_close n b v w ( by aesop ) ; aesop
  refine le_trans ( Finset.card_le_card h_partition ) ?_
  exact le_trans ( Finset.card_union_le _ _ )
    ( add_le_add
      ( le_trans ( Finset.card_union_le _ _ )
        ( add_le_add ( card_neighbors_same_time n b v )
                     ( card_neighbors_next_time n b v ) ) )
      ( card_neighbors_prev_time n b v ) )
