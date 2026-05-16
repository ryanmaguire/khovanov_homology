/-
  Fixed-Strand Conjecture for 4-strand braids: Small cases (length ≤ 2)
  
  Strategy: Define a BFS-based Bool connectivity check on the complement graph,
  prove it implies is_connected for the independence complex,
  then use native_decide on the Bool check.
-/
import FixedStrand.Defs

set_option maxHeartbeats 800000

/-! ## BFS connectivity check -/

/-- All vertices of the 4-strand augmented rhomboid graph. -/
def allVertices4 (b : Braid 4) : List (Fin 4 × Fin (b.length + 1)) :=
  (List.finRange 4).flatMap fun s => (List.finRange (b.length + 1)).map fun t => (s, t)

/-- Bool adjacency for the augmented rhomboid graph. -/
def augAdjBool (b : Braid 4) (u v : Fin 4 × Fin (b.length + 1)) : Bool :=
  (u.1 == v.1 && (u.2.val + 1 == v.2.val || v.2.val + 1 == u.2.val)) ||
  (u.2 == v.2 && (u.1.val + 1 == v.1.val || v.1.val + 1 == u.1.val)) ||
  (List.finRange b.length).any fun (t : Fin b.length) =>
    let i := (b.get t).val
    (u.1.val == i && u.2.val == t.val && v.1.val == i + 1 && v.2.val == t.val + 1) ||
    (u.1.val == i + 1 && u.2.val == t.val && v.1.val == i && v.2.val == t.val + 1) ||
    (v.1.val == i && v.2.val == t.val && u.1.val == i + 1 && u.2.val == t.val + 1) ||
    (v.1.val == i + 1 && v.2.val == t.val && u.1.val == i && u.2.val == t.val + 1)

/-- Bool complement adjacency: u ≠ v and not G-adjacent. -/
def complAdjBool (b : Braid 4) (u v : Fin 4 × Fin (b.length + 1)) : Bool :=
  u != v && !augAdjBool b u v

/-- BFS connectivity check on the complement graph.
    Returns true iff all vertices are reachable from the first vertex. -/
def isComplConnected (b : Braid 4) : Bool :=
  let verts := allVertices4 b
  match verts with
  | [] => true
  | start :: _ =>
    let visited := (List.range verts.length).foldl (fun vis _ =>
      let newVerts := verts.filter fun v => !vis.contains v && vis.any fun u => complAdjBool b u v
      (vis ++ newVerts).eraseDups
    ) [start]
    verts.all fun v => visited.contains v

/-! ## Correctness of BFS check -/

/-
augAdjBool captures the symmetrized adjacency relation (both directions).
-/
lemma augAdjBool_iff_sym (b : Braid 4) (u v : Fin 4 × Fin (b.length + 1)) :
    augAdjBool b u v = true ↔ (AugRhomboidAdj b u v ∨ AugRhomboidAdj b v u) := by
  constructor <;> intro h <;> simp_all +decide [ AugRhomboidAdj, augAdjBool ];
  · grind;
  · grind

/-
complAdjBool correctly characterizes complement adjacency.
-/
lemma complAdjBool_iff (b : Braid 4) (u v : Fin 4 × Fin (b.length + 1)) :
    complAdjBool b u v = true ↔ (u ≠ v ∧ ¬(AugmentedRhomboidGraph b).Adj u v) := by
  unfold complAdjBool AugmentedRhomboidGraph;
  have := augAdjBool_iff_sym b u v; aesop;

/-
If isComplConnected returns true, then the independence complex is connected.
    The key idea: BFS finds that every vertex is reachable from vertex 0 in the complement.
    A complement path corresponds to a path in the independence complex.
-/
lemma isComplConnected_implies_indep_connected (b : Braid 4)
    (h : isComplConnected b = true) :
    is_connected (IndependenceComplex (AugmentedRhomboidGraph b)) := by
  -- By definition of `isComplConnected`, if `isComplConnected b` is true, then the complement graph is connected.
  have h_compl_connected : SimpleGraph.Connected (SimpleGraph.comap (fun v => v) (AugmentedRhomboidGraph b)ᶜ) := by
    rw [ SimpleGraph.connected_iff_exists_forall_reachable ];
    use (0, 0);
    intro w
    have h_reachable : w ∈ (List.range (allVertices4 b).length).foldl (fun (vis : List (Fin 4 × Fin (b.length + 1))) (x : ℕ) =>
      let newVerts := (allVertices4 b).filter (fun (v : Fin 4 × Fin (b.length + 1)) => !vis.contains v && vis.any (fun (u : Fin 4 × Fin (b.length + 1)) => complAdjBool b u v))
      (vis ++ newVerts).eraseDups
    ) [(0, 0)] := by
      unfold isComplConnected at h;
      unfold allVertices4 at *; simp_all +decide [ List.range_succ_eq_map ] ;
      simp_all +decide [ List.finRange ];
      rcases w with ⟨ i, j ⟩ ; fin_cases i <;> simp_all +decide ;
      · exact if hj : j = 0 then hj.symm ▸ h.1.1 else h.1.2 j hj;
      · grind;
      · grind +locals;
      · grind;
    -- By definition of `foldl`, we can show that every vertex in the list is reachable from `(0, 0)`.
    have h_reachable : ∀ (vis : List (Fin 4 × Fin (b.length + 1))) (x : ℕ), (∀ v ∈ vis, SimpleGraph.Reachable (SimpleGraph.comap (fun v => v) (AugmentedRhomboidGraph b)ᶜ) (0, 0) v) → (∀ v ∈ List.filter (fun v => !vis.contains v && vis.any fun u => complAdjBool b u v) (allVertices4 b), SimpleGraph.Reachable (SimpleGraph.comap (fun v => v) (AugmentedRhomboidGraph b)ᶜ) (0, 0) v) := by
      intros vis x hvis v hv
      obtain ⟨u, hu⟩ : ∃ u ∈ vis, complAdjBool b u v := by
        grind;
      have h_reachable : SimpleGraph.Adj (SimpleGraph.comap (fun v => v) (AugmentedRhomboidGraph b)ᶜ) u v := by
        simp_all +decide [ complAdjBool_iff ];
      exact SimpleGraph.Reachable.trans ( hvis u hu.1 ) ( SimpleGraph.Adj.reachable h_reachable );
    have h_reachable : ∀ (vis : List (Fin 4 × Fin (b.length + 1))) (x : ℕ), (∀ v ∈ vis, SimpleGraph.Reachable (SimpleGraph.comap (fun v => v) (AugmentedRhomboidGraph b)ᶜ) (0, 0) v) → (∀ v ∈ List.eraseDups (vis ++ List.filter (fun v => !vis.contains v && vis.any fun u => complAdjBool b u v) (allVertices4 b)), SimpleGraph.Reachable (SimpleGraph.comap (fun v => v) (AugmentedRhomboidGraph b)ᶜ) (0, 0) v) := by
      intros vis x hvis v hv
      have h_eraseDups : ∀ (l : List (Fin 4 × Fin (b.length + 1))), (∀ v ∈ l, SimpleGraph.Reachable (SimpleGraph.comap (fun v => v) (AugmentedRhomboidGraph b)ᶜ) (0, 0) v) → (∀ v ∈ List.eraseDups l, SimpleGraph.Reachable (SimpleGraph.comap (fun v => v) (AugmentedRhomboidGraph b)ᶜ) (0, 0) v) := by
        intros l hl v hv;
        induction' l using List.reverseRecOn with l ih;
        · contradiction;
        · simp_all +decide [ List.eraseDups_append ];
          simp_all +decide [ List.removeAll ];
          grind;
      exact h_eraseDups _ ( fun v hv => by cases List.mem_append.mp hv <;> [ exact hvis _ ‹_›; exact h_reachable _ x hvis _ ‹_› ] ) _ hv;
    have h_reachable : ∀ (l : List ℕ), (∀ v ∈ List.foldl (fun (vis : List (Fin 4 × Fin (b.length + 1))) (x : ℕ) =>
      let newVerts := List.filter (fun v => !vis.contains v && vis.any fun u => complAdjBool b u v) (allVertices4 b)
      (vis ++ newVerts).eraseDups
    ) [(0, 0)] l, SimpleGraph.Reachable (SimpleGraph.comap (fun v => v) (AugmentedRhomboidGraph b)ᶜ) (0, 0) v) := by
      intro l
      induction' l using List.reverseRecOn with l ih;
      · simp +decide [ SimpleGraph.Reachable.refl ];
      · convert h_reachable _ _ ‹_› using 1;
        · simp +decide [ List.foldl_append ];
        · exact 0;
    exact h_reachable _ _ ‹_›;
  intro u v hu hv;
  obtain ⟨ p, hp ⟩ := h_compl_connected u v;
  · use [u]; simp [hu];
    exact List.isChain_singleton _;
  · rename_i w hw₁ hw₂;
    have h_path : ∀ {x y : Fin 4 × Fin (List.length b + 1)}, (SimpleGraph.comap (fun v => v) (AugmentedRhomboidGraph b)ᶜ).Adj x y → {x, y} ∈ (IndependenceComplex (AugmentedRhomboidGraph b)).faces := by
      simp +decide [ IndependenceComplex ];
      rintro a b c d h₁ h₂ e f ( ⟨ rfl, rfl ⟩ | ⟨ rfl, rfl ⟩ ) g h ( ⟨ rfl, rfl ⟩ | ⟨ rfl, rfl ⟩ ) <;> simp_all +decide [ SimpleGraph.adj_comm ];
    have h_path : ∀ {x y : Fin 4 × Fin (List.length b + 1)}, (SimpleGraph.comap (fun v => v) (AugmentedRhomboidGraph b)ᶜ).Walk x y → ∃ p : List (Fin 4 × Fin (List.length b + 1)), p ≠ [] ∧ p.head? = some x ∧ p.getLast? = some y ∧ List.Chain' (fun a b_1 => {a, b_1} ∈ (IndependenceComplex (AugmentedRhomboidGraph b)).faces) p := by
      intro x y p;
      induction' p with x y p ih;
      · use [x]; simp;
        exact List.isChain_singleton _;
      · rename_i h₁ h₂ h₃;
        obtain ⟨ p, hp₁, hp₂, hp₃, hp₄ ⟩ := h₃; use y :: p; simp_all +decide [ List.chain'_cons' ] ;
        exact ⟨ by cases p <;> aesop, List.chain'_cons'.mpr ⟨ by aesop, hp₄ ⟩ ⟩;
    exact h_path ( SimpleGraph.Walk.cons hw₁ hw₂ )

/-! ## Small braid enumeration -/

/-- All 4-strand braids of length ≤ 2. -/
def smallBraids4 : List (Braid 4) :=
  [[]] ++
  [[0], [1], [2]] ++
  [[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]]

/-
Every braid of length < 3 is in smallBraids4.
-/
lemma mem_smallBraids4_of_length_lt_3 (b : Braid 4) (h : b.length < 3) :
    b ∈ smallBraids4 := by
  interval_cases _ : b.length <;> simp_all +decide [ smallBraids4 ];
  · rcases b with ( _ | ⟨ a, _ | ⟨ b, _ | ⟨ c, _ | b ⟩ ⟩ ⟩ ) <;> simp_all +decide;
    decide +revert;
  · rcases b with ( _ | ⟨ a, _ | ⟨ b, _ | b ⟩ ⟩ ) <;> simp_all +decide;
    native_decide +revert

/-- All small braids pass the BFS connectivity check. -/
lemma allSmallBraidsConnectedBool :
    smallBraids4.all (fun b => isComplConnected b) = true := by
  native_decide

/-- All 4-strand braids of length ≤ 2 have connected independence complex. -/
theorem small_braids_connected (b : Braid 4) (h : b.length < 3) :
    is_connected (IndependenceComplex (AugmentedRhomboidGraph b)) := by
  have hmem := mem_smallBraids4_of_length_lt_3 b h
  have hall := allSmallBraidsConnectedBool
  have : isComplConnected b = true := by
    rw [List.all_eq_true] at hall
    exact hall b hmem
  exact isComplConnected_implies_indep_connected b this
