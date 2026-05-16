/-
  Fixed-Strand Conjecture for 4-strand braids: Large cases (length ≥ 3)
  
  Strategy: Show that each vertex in the augmented rhomboid graph has degree ≤ 6.
  For length ≥ 3, the vertex set has size 4(L+1) ≥ 16, so any two vertices have 
  a common non-neighbor (diameter ≤ 2 in the complement), which implies connectivity.
-/
import FixedStrand.Defs
import FixedStrand.Connectivity

set_option maxHeartbeats 800000

open scoped Classical

variable {V : Type*} [DecidableEq V]

/-! ## Degree bound -/

/-
The degree of a vertex in the augmented rhomboid graph is at most 6.
-/
lemma aug_degree_le_6 (b : Braid 4) (v : Fin 4 × Fin (b.length + 1)) :
    (Finset.univ.filter fun w => (AugmentedRhomboidGraph b).Adj v w).card ≤ 6 := by
  refine' le_trans ( Finset.card_le_card _ ) _;
  exact { ( v.1, v.2 + 1 ), ( v.1, v.2 - 1 ), ( v.1 + 1, v.2 ), ( v.1 - 1, v.2 ) } ∪ Finset.image ( fun t : Fin b.length => if v.1 = ( b.get t ).val then ( ( v.1 + 1, t.succ ) ) else if v.1 = ( b.get t ).val + 1 then ( ( v.1 - 1, t.succ ) ) else ( v.1, v.2 ) ) ( Finset.univ.filter fun t : Fin b.length => v.2 = t.castSucc ) ∪ Finset.image ( fun t : Fin b.length => if v.1 = ( b.get t ).val + 1 then ( ( v.1 - 1, t.castSucc ) ) else if v.1 = ( b.get t ).val then ( ( v.1 + 1, t.castSucc ) ) else ( v.1, v.2 ) ) ( Finset.univ.filter fun t : Fin b.length => v.2 = t.succ );
  · intro w hw; simp_all +decide [ Finset.subset_iff ] ;
    rcases hw with ⟨ hw₁, hw₂ ⟩ ; rcases hw₂ with ( ( hw₂ | hw₂ | hw₂ ) | hw₂ | hw₂ | hw₂ ) <;> simp_all +decide [ Fin.ext_iff, Prod.ext_iff ] ;
    all_goals norm_num [ Fin.val_add, Fin.val_sub ] at *;
    any_goals omega;
    · cases hw₂.2 <;> simp_all +decide [ Nat.mod_eq_of_lt ];
      · exact Or.inl ( by rw [ Nat.mod_eq_of_lt ] ; linarith [ Fin.is_lt v.2, Fin.is_lt w.2 ] );
      · refine' Or.inr <| Or.inl <| _;
        rcases b with ⟨ _ | _ | b ⟩ <;> simp_all +arith +decide [ Nat.mod_eq_of_lt ];
        rw [ Nat.mod_eq_sub_mod ] <;> norm_num [ ← ‹ ( w.2 : ℕ ) + 1 = v.2 › ];
        · rw [ Nat.mod_eq_of_lt ] <;> omega;
        · linarith;
    · rcases hw₂ with ⟨ t, ht | ht ⟩ <;> simp_all +decide [ Fin.ext_iff ];
      · grind;
      · grind +revert;
    · rcases hw₂.2 with ( hw₂ | hw₂ ) <;> simp_all +decide [ Nat.mod_eq_of_lt ];
      · refine' Or.inr <| Or.inl <| _;
        rcases b with ( _ | ⟨ _, _ | b ⟩ ) <;> simp_all +arith +decide [ Nat.mod_eq_of_lt ];
        · native_decide +revert;
        · rw [ Nat.mod_eq_sub_mod ] <;> norm_num [ hw₂.symm ];
          · rw [ Nat.mod_eq_of_lt ] <;> omega;
          · linarith;
      · exact Or.inl ( by rw [ Nat.mod_eq_of_lt ] ; linarith [ Fin.is_lt v.2, Fin.is_lt w.2 ] );
    · rcases hw₂ with ⟨ t, ht | ht ⟩ <;> simp_all +decide [ Fin.ext_iff ];
      · grind +locals;
      · grind;
  · refine' le_trans ( Finset.card_union_le _ _ ) _;
    refine' le_trans ( add_le_add ( Finset.card_union_le _ _ ) ( Finset.card_image_le ) ) _;
    refine' le_trans ( add_le_add_three ( Finset.card_insert_le _ _ ) ( Finset.card_image_le ) ( Finset.card_le_one.mpr _ ) ) _;
    · aesop;
    · refine' le_trans ( add_le_add ( add_le_add ( add_le_add ( Finset.card_insert_le _ _ ) le_rfl ) ( Finset.card_le_one.mpr _ ) ) le_rfl ) _ <;> norm_num;
      · aesop;
      · grind

/-! ## Complement connectivity from degree bound -/

/-
If every vertex has degree ≤ k in G, and |V| > 2k + 2, then Gᶜ is connected
    (any two vertices share a common non-neighbor).
-/
lemma compl_connected_of_low_degree {W : Type*} [DecidableEq W] [Fintype W]
    (G : SimpleGraph W) [DecidableRel G.Adj]
    (k : ℕ) (hk : ∀ v, (Finset.univ.filter fun w => G.Adj v w).card ≤ k)
    (hV : 2 * k + 2 < Fintype.card W) :
    Gᶜ.Connected := by
  have h_g_compl_connected : ∀ u v : W, u ≠ v → ¬G.Adj u v ∨ ∃ w : W, w ≠ u ∧ w ≠ v ∧ ¬G.Adj u w ∧ ¬G.Adj v w := by
    intro u v huv
    by_contra h_contra
    push_neg at h_contra;
    have h_card : (Finset.univ \ ({u, v} : Finset W)).card ≤ (Finset.filter (fun w => G.Adj u w) Finset.univ).card + (Finset.filter (fun w => G.Adj v w) Finset.univ).card := by
      have h_card : (Finset.univ \ ({u, v} : Finset W)) ⊆ (Finset.filter (fun w => G.Adj u w) Finset.univ) ∪ (Finset.filter (fun w => G.Adj v w) Finset.univ) := by
        grind;
      exact le_trans ( Finset.card_le_card h_card ) ( Finset.card_union_le _ _ );
    simp_all +decide [ Finset.card_sdiff, Finset.card_singleton ];
    linarith [ hk u, hk v ];
  have h_g_compl_connected : ∀ u v : W, u ≠ v → Gᶜ.Reachable u v := by
    intro u v huv
    obtain h | ⟨w, hw₁, hw₂, hw₃, hw₄⟩ := h_g_compl_connected u v huv;
    · exact SimpleGraph.Adj.reachable ( by aesop );
    · have h_g_compl_connected : Gᶜ.Reachable u w ∧ Gᶜ.Reachable w v := by
        exact ⟨ SimpleGraph.Adj.reachable ( by tauto ), SimpleGraph.Adj.reachable ( by tauto ) ⟩;
      exact h_g_compl_connected.1.trans h_g_compl_connected.2;
  rcases isEmpty_or_nonempty W with h | h;
  · aesop;
  · exact SimpleGraph.Connected.mk fun u v => if huv : u = v then huv.symm ▸ SimpleGraph.Reachable.refl _ else h_g_compl_connected u v huv

/-
For 4-strand braids of length ≥ 3, the independence complex is connected.
-/
theorem large_braids_connected (b : Braid 4) (h : b.length ≥ 3) :
    is_connected (IndependenceComplex (AugmentedRhomboidGraph b)) := by
  apply indep_connected_of_compl_connected;
  apply compl_connected_of_low_degree;
  exact fun v => aug_degree_le_6 b v;
  norm_num; linarith
