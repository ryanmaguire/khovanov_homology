/-
  Fixed-Strand Conjecture for 4-strand braids: Collapse machinery
-/
import FixedStrand.Defs

set_option maxHeartbeats 800000

open Classical

variable {V : Type*} [DecidableEq V]

/-! ## Elementary collapse preserves connectivity -/

/-- If a face is in the elementary collapse, it is in the original complex. -/
lemma elementary_collapse_faces_subset (K : AbstractSimplicialComplex V) (σ τ : Finset V)
    (h : is_free_pair K σ τ) :
    (elementary_collapse K σ τ h).faces ⊆ K.faces := by
  intro f hf; exact hf.1

/-- Case σ.card ≥ 2: edges can be rerouted through τ's extra vertex. -/
lemma elementary_collapse_preserves_connected_high_dim
    (K : AbstractSimplicialComplex V) (σ τ : Finset V)
    (h : is_free_pair K σ τ) (hconn : is_connected K) (hcard : σ.card ≥ 2) :
    is_connected (elementary_collapse K σ τ h) := by
  have edge_bypass : ∀ a b, {a, b} ∈ K.faces → ¬ ({a, b} ∈ (elementary_collapse K σ τ h).faces) → ∃ c, {a, c} ∈ (elementary_collapse K σ τ h).faces ∧ {c, b} ∈ (elementary_collapse K σ τ h).faces := by
    intro a b hab hnot
    have h_eq : {a, b} = σ := by
      contrapose! hnot; simp_all +decide [ elementary_collapse ] ;
      intro h; have := h ▸ h; simp_all +decide [ is_free_pair ] ;
      grind;
    have hτ_card : τ.card ≥ 3 := by linarith [ h.2.2.2.1 ];
    obtain ⟨c, hc⟩ : ∃ c ∈ τ, c ∉ σ := by
      exact Finset.not_subset.mp fun h' => by have := Finset.card_le_card h'; linarith [ h.2.2.2.1 ] ;
    have h_ac_bc : {a, c} ∈ K.faces ∧ {c, b} ∈ K.faces := by
      have h_ac_bc : {a, c} ⊆ τ ∧ {c, b} ⊆ τ := by
        have := h.2.1; simp_all +decide [ Finset.subset_iff ] ;
        have := h.2.2.1; aesop;
      exact ⟨ K.down_closed _ h.2.1 _ h_ac_bc.1, K.down_closed _ h.2.1 _ h_ac_bc.2 ⟩;
    refine' ⟨ c, _, _ ⟩ <;> simp_all +decide [ elementary_collapse ];
    · grind;
    · grind;
  intro u v hu hv;
  obtain ⟨ p, hp ⟩ := hconn u v ( elementary_collapse_faces_subset K σ τ h hu ) ( elementary_collapse_faces_subset K σ τ h hv );
  have h_chain_bypass : ∀ (l : List V), l ≠ [] → List.Chain' (fun a b => {a, b} ∈ K.faces) l → ∃ l' : List V, l'.head? = l.head? ∧ l'.getLast? = l.getLast? ∧ List.Chain' (fun a b => {a, b} ∈ (elementary_collapse K σ τ h).faces) l' := by
    intro l hl hl_chain
    induction' l with a l ih;
    · contradiction;
    · rcases l with ( _ | ⟨ b, l ⟩ ) <;> simp_all +decide [ List.chain'_cons' ];
      · exact ⟨ [ a ], rfl, rfl, List.isChain_singleton _ ⟩;
      · by_cases hab : {a, b} ∈ (elementary_collapse K σ τ h).faces;
        · obtain ⟨ l', hl' ⟩ := ih ( List.chain'_cons'.1 hl_chain |>.2 );
          use a :: l';
          cases l' <;> simp_all +decide [ List.chain'_cons' ];
          exact ⟨ by aesop, List.isChain_cons_cons.mpr ⟨ hab, by aesop ⟩ ⟩;
        · obtain ⟨ c, hc₁, hc₂ ⟩ := edge_bypass a b ( by
            exact List.isChain_cons_cons.mp hl_chain |>.1 ) hab;
          obtain ⟨ l', hl'₁, hl'₂, hl'₃ ⟩ := ih ( List.chain'_cons'.1 hl_chain |>.2 );
          use a :: c :: l';
          cases l' <;> simp_all +decide [ List.chain'_cons' ];
          exact List.isChain_cons_cons.mpr ⟨ hc₁, List.isChain_cons_cons.mpr ⟨ hc₂, hl'₃ ⟩ ⟩;
  grind

/-
Case σ.card = 1: removing a pendant vertex.
-/
lemma elementary_collapse_preserves_connected_dim1
    (K : AbstractSimplicialComplex V) (σ τ : Finset V)
    (h : is_free_pair K σ τ) (hconn : is_connected K) (hcard : σ.card = 1) :
    is_connected (elementary_collapse K σ τ h) := by
  intro u v hu hv;
  obtain ⟨w, hw⟩ : ∃ w, σ = {w} := by
    exact Finset.card_eq_one.mp hcard;
  -- Since $u$ and $v$ are in the original complex and not equal to $w$, there exists a path $p$ in $K$ from $u$ to $v$ that avoids $w$.
  obtain ⟨p, hp⟩ : ∃ p : List V, p ≠ [] ∧ p.head? = some u ∧ p.getLast? = some v ∧ List.Chain' (fun a b => {a, b} ∈ K.faces) p ∧ ∀ x ∈ p, x ≠ w := by
    obtain ⟨ p, hp ⟩ := hconn u v ( by
      exact elementary_collapse_faces_subset K σ τ h hu ) ( by
      exact elementary_collapse_faces_subset K σ τ h hv );
    -- Replace any occurrence of $w$ in $p$ with the unique vertex $w'$ adjacent to $w$ in $\tau$.
    obtain ⟨w', hw'⟩ : ∃ w', τ = {w, w'} := by
      have hτ_card : τ.card = 2 := by
        linarith [ h.2.2.2.1 ];
      have hτ_subset : {w} ⊆ τ := by
        have := h.2.2.1; aesop;
      rw [ Finset.card_eq_two ] at hτ_card; obtain ⟨ x, y, hxy ⟩ := hτ_card; use if x = w then y else x; aesop;
    -- Replace any occurrence of $w$ in $p$ with $w'$.
    use p.map (fun x => if x = w then w' else x);
    refine' ⟨ _, _, _, _, _ ⟩;
    · aesop;
    · unfold elementary_collapse at hu; aesop;
    · unfold elementary_collapse at hv; aesop;
    · have h_chain : ∀ a b : V, {a, b} ∈ K.faces → {if a = w then w' else a, if b = w then w' else b} ∈ K.faces := by
        intro a b hab
        by_cases ha : a = w
        by_cases hb : b = w;
        · have := h.2.1; simp_all +decide [ Finset.subset_iff ] ;
          exact K.down_closed _ this _ ( by simp +decide );
        · have := h.2.2.2.2;
          grind +suggestions;
        · by_cases hb : b = w <;> simp_all +decide;
          have := h.2.2.2.2 { a, w } hab; simp_all +decide [ Finset.subset_iff ] ;
          grind +suggestions;
      have h_chain : ∀ {l : List V}, List.Chain' (fun a b => {a, b} ∈ K.faces) l → List.Chain' (fun a b => {a, b} ∈ K.faces) (List.map (fun x => if x = w then w' else x) l) := by
        intro l hl; induction hl <;> simp_all +decide [ List.isChain_cons_cons ] ;
        · exact List.isChain_nil;
        · exact List.isChain_singleton _;
        · exact List.IsChain.cons_cons ( h_chain _ _ ‹_› ) ‹_›;
      exact h_chain hp.2.2.2;
    · simp +contextual [ hw' ];
      intro a ha; split_ifs <;> simp_all +decide ;
      have := h.2.2.2.1; simp_all +decide ;
      rintro rfl; simp +decide at this;
  refine' ⟨ p, hp.1, hp.2.1, hp.2.2.1, _ ⟩;
  have h_edges : ∀ a b, {a, b} ∈ K.faces → a ≠ w → b ≠ w → {a, b} ∈ (elementary_collapse K σ τ h).faces := by
    intro a b hab ha hb
    simp [elementary_collapse, hw];
    refine' ⟨ hab, _, _ ⟩ <;> intro H <;> simp_all +decide [ Finset.eq_singleton_iff_unique_mem ];
    have := h.2.2.1; simp_all +decide [ Finset.subset_iff ] ;
    grind;
  have := List.isChain_iff_get.mp hp.2.2.2.1;
  refine' List.isChain_iff_get.mpr _;
  exact fun i => h_edges _ _ ( this i ) ( hp.2.2.2.2 _ ( by simp ) ) ( hp.2.2.2.2 _ ( by simp ) )

/-- Case σ.card = 0: complex becomes empty (vacuous). -/
lemma elementary_collapse_preserves_connected_dim0
    (K : AbstractSimplicialComplex V) (σ τ : Finset V)
    (h : is_free_pair K σ τ) (hconn : is_connected K) (hcard : σ.card = 0) :
    is_connected (elementary_collapse K σ τ h) := by
  have h_faces_subset : K.faces ⊆ {∅, τ} := by cases h; aesop
  intro u v hu hv
  have hu_K := hu.1
  have hv_K := hv.1
  have hu_set : ({u} : Finset V) ∈ ({∅, τ} : Set (Finset V)) := h_faces_subset hu_K
  simp at hu_set
  -- {u} = ∅ is impossible, so {u} = τ, meaning τ = {u}. But then τ.card = 1.
  -- Since τ.card = σ.card + 1 = 1, this is fine.
  -- But {u} should not be in L.faces (since {u} = τ).
  -- Since hu says {u} ∈ L.faces = K.faces \ {σ, τ}, we need {u} ≠ σ and {u} ≠ τ.
  -- But {u} = τ, contradiction.
  exfalso
  rw [hu_set] at hu
  simp only [elementary_collapse, Set.mem_diff, Set.mem_insert_iff, Set.mem_singleton_iff] at hu
  tauto

/-- Elementary collapse preserves connectivity. -/
lemma elementary_collapse_preserves_connected (K : AbstractSimplicialComplex V) (σ τ : Finset V)
    (h : is_free_pair K σ τ) (hconn : is_connected K) :
    is_connected (elementary_collapse K σ τ h) := by
  by_cases h2 : σ.card ≥ 2
  · exact elementary_collapse_preserves_connected_high_dim K σ τ h hconn h2
  · by_cases h1 : σ.card = 1
    · exact elementary_collapse_preserves_connected_dim1 K σ τ h hconn h1
    · exact elementary_collapse_preserves_connected_dim0 K σ τ h hconn (by omega)

/-- Collapses preserve connectivity. -/
lemma Collapses_preserves_connected {K L : AbstractSimplicialComplex V}
    (hcol : Collapses K L) (hconn : is_connected K) : is_connected L := by
  induction hcol with
  | refl => exact hconn
  | step σ τ h _ ih =>
    exact elementary_collapse_preserves_connected _ σ τ h ih

/-! ## collapse_all -/

noncomputable def collapse_all_aux (K : AbstractSimplicialComplex V) : ℕ → AbstractSimplicialComplex V
  | 0 => K
  | fuel + 1 =>
    if h : ∃ σ τ, is_free_pair K σ τ then
      collapse_all_aux (elementary_collapse K h.choose h.choose_spec.choose h.choose_spec.choose_spec) fuel
    else K

noncomputable def collapse_all [Fintype V] (K : AbstractSimplicialComplex V) : AbstractSimplicialComplex V :=
  collapse_all_aux K (2 ^ Fintype.card (Finset V))

lemma collapse_all_aux_collapses (K : AbstractSimplicialComplex V) (fuel : ℕ) :
    Collapses K (collapse_all_aux K fuel) := by
  induction fuel generalizing K with
  | zero => exact Collapses.refl K
  | succ n ih =>
    simp only [collapse_all_aux]
    split
    · next h =>
      exact Collapses_trans
        (Collapses.step _ _ h.choose_spec.choose_spec (Collapses.refl K))
        (ih _)
    · exact Collapses.refl K

lemma collapse_all_collapses [Fintype V] (K : AbstractSimplicialComplex V) :
    Collapses K (collapse_all K) :=
  collapse_all_aux_collapses K _

lemma elementary_collapse_card_lt (K : AbstractSimplicialComplex V) (σ τ : Finset V)
    (h : is_free_pair K σ τ) (hfin : K.faces.Finite) :
    (hfin.subset (elementary_collapse_faces_subset K σ τ h)).toFinset.card <
      hfin.toFinset.card := by
  have h_card : (hfin.toFinset \ {σ, τ}).card < hfin.toFinset.card := by
    refine' Finset.card_lt_card _;
    simp_all +decide [ Finset.ssubset_def, Finset.subset_iff ];
    exact ⟨ σ, h.1, by tauto ⟩;
  convert h_card using 2;
  ext; simp [elementary_collapse]

lemma collapse_all_aux_no_free_pairs (K : AbstractSimplicialComplex V) (fuel : ℕ)
    (hfin : K.faces.Finite) (hfuel : hfin.toFinset.card ≤ fuel) :
    has_no_free_pairs (collapse_all_aux K fuel) := by
  induction' fuel with fuel ih generalizing K;
  · simp_all +decide [ Finset.ext_iff, has_no_free_pairs ];
    exact fun σ τ h => hfuel σ h.1;
  · by_cases h : ∃ σ τ, is_free_pair K σ τ;
    · convert ih ( elementary_collapse K h.choose h.choose_spec.choose h.choose_spec.choose_spec ) _ _ using 1;
      rotate_left;
      exact Set.Finite.subset hfin ( elementary_collapse_faces_subset K _ _ h.choose_spec.choose_spec );
      · exact Nat.le_of_lt_succ ( lt_of_lt_of_le ( elementary_collapse_card_lt K _ _ h.choose_spec.choose_spec hfin ) hfuel );
      · exact dif_pos h;
    · unfold collapse_all_aux; aesop;

lemma collapse_all_has_no_free_pairs [Fintype V] (K : AbstractSimplicialComplex V)
    (hfin : K.faces.Finite) : has_no_free_pairs (collapse_all K) := by
  apply collapse_all_aux_no_free_pairs;
  refine' le_trans ( Finset.card_le_univ _ ) _;
  · exact le_of_lt ( Nat.recOn ( Fintype.card ( Finset V ) ) ( by decide ) fun n ihn => by rw [ pow_succ' ] ; linarith );
  · exact hfin
