/-
  Fixed-Strand Conjecture for 4-strand braids: Connectivity
  
  Key result: is_connected (IndependenceComplex G) ↔ Gᶜ.Connected (for nonempty finite types)
-/
import FixedStrand.Defs

set_option maxHeartbeats 800000

open scoped Classical

variable {V : Type*} [DecidableEq V] [Fintype V]

/-! ## Independence complex connectivity ↔ complement graph connectivity -/

/-- If two vertices are non-adjacent in G, then {u, v} is a face of the independence complex. -/
lemma indep_edge_of_not_adj (G : SimpleGraph V) (u v : V) (huv : u ≠ v) (h : ¬G.Adj u v) :
    ({u, v} : Finset V) ∈ (IndependenceComplex G).faces := by
  simp only [IndependenceComplex, Set.mem_setOf_eq]
  intro a ha b hb hab
  simp only [Finset.mem_insert, Finset.mem_singleton] at ha hb
  rcases ha with rfl | rfl <;> rcases hb with rfl | rfl
  · exact absurd rfl hab
  · exact h
  · exact fun hadj => h (G.symm hadj)
  · exact absurd rfl hab

/-- Any singleton is a face of the independence complex. -/
lemma indep_singleton (G : SimpleGraph V) (v : V) :
    ({v} : Finset V) ∈ (IndependenceComplex G).faces := by
  simp only [IndependenceComplex, Set.mem_setOf_eq]
  intro a ha b hb hab
  simp only [Finset.mem_singleton] at ha hb
  subst ha; subst hb; exact absurd rfl hab

/-
If Gᶜ is connected, then the independence complex of G is connected.
-/
lemma indep_connected_of_compl_connected (G : SimpleGraph V) [DecidableRel G.Adj]
    (hconn : Gᶜ.Connected) : is_connected (IndependenceComplex G) := by
  intro u v;
  -- Since $Gᶜ$ is connected, there exists a walk from $u$ to $v$ in $Gᶜ$.
  obtain ⟨w, hw⟩ : ∃ w : Gᶜ.Walk u v, True := by
    have := hconn u v; aesop;
  generalize_proofs at *; (
  intro hu hv
  use w.support
  generalize_proofs at *; (
  refine' ⟨ _, _, _, _ ⟩ <;> induction w <;> simp_all +decide [ List.chain'_cons' ];
  · cases ‹Gᶜ.Walk _ _› <;> simp_all +decide [ List.getLast? ];
  · exact List.isChain_singleton _;
  · rename_i u v w hw ih; induction hw <;> simp_all +decide [ List.chain'_cons' ] ;
    · exact List.isChain_cons_cons.mpr ⟨ by simpa [ SimpleGraph.adj_comm ] using indep_edge_of_not_adj G _ _ w.1 w.2, ih ⟩;
    · exact List.isChain_cons_cons.mpr ⟨ by simpa [ Finset.pair_comm ] using indep_edge_of_not_adj G _ _ ( by tauto ) ( by tauto ), ih ( by simpa [ Finset.pair_comm ] using indep_singleton G _ ) ⟩))

/-
If the independence complex of G is connected, then Gᶜ is connected.
-/
lemma compl_connected_of_indep_connected (G : SimpleGraph V) [DecidableRel G.Adj]
    (hconn : is_connected (IndependenceComplex G)) [Nonempty V] : Gᶜ.Connected := by
  refine' ⟨ fun u v => _ ⟩;
  obtain ⟨ p, hp ⟩ := hconn u v ( indep_singleton G u ) ( indep_singleton G v );
  induction' p with a p ih generalizing u v <;> simp_all +decide [ List.chain'_cons' ];
  rcases p with ( _ | ⟨ b, p ⟩ ) <;> simp_all +decide [ List.getLast? ];
  · aesop;
  · have h_adj : ¬G.Adj a b := by
      have := List.isChain_cons_cons.mp hp.2.2; simp_all +decide [ IndependenceComplex ] ;
      by_cases h : u = b <;> simp_all +decide [ SimpleGraph.adj_comm ];
    have h_adj : Gᶜ.Reachable a b := by
      by_cases hab : a = b <;> simp_all +decide [ SimpleGraph.Reachable ];
      · exact ⟨ SimpleGraph.Walk.nil ⟩;
      · exact ⟨ SimpleGraph.Walk.cons ( by aesop ) SimpleGraph.Walk.nil ⟩;
    convert h_adj.trans ( ih b rfl ( List.chain'_cons'.1 hp.2.2 |>.2 ) ) using 1;
    · exact hp.1.symm;
    · exact hp.2.1.symm

/-- For nonempty finite types, independence complex connectivity ↔ complement connectivity. -/
lemma indep_connected_iff_compl_connected (G : SimpleGraph V) [DecidableRel G.Adj] [Nonempty V] :
    is_connected (IndependenceComplex G) ↔ Gᶜ.Connected :=
  ⟨fun h => compl_connected_of_indep_connected G h, fun h => indep_connected_of_compl_connected G h⟩
