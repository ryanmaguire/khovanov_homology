/-
  Fixed-Strand Conjecture for 4-strand braids: Definitions
-/
import Mathlib

set_option maxHeartbeats 400000

/-! ## Basic Definitions -/

/-- A braid on n strands is a list of crossings, each a value in Fin (n-1). -/
def Braid (n : ℕ) := List (Fin (n - 1))

/-- The adjacency relation for the augmented rhomboid graph. -/
def AugRhomboidAdj {n : ℕ} (b : Braid n) (u v : Fin n × Fin (b.length + 1)) : Prop :=
  -- Horizontal edges (same strand, adjacent time steps)
  (u.1 = v.1 ∧ (u.2.val + 1 = v.2.val ∨ v.2.val + 1 = u.2.val)) ∨
  -- Vertical edges (same time step, adjacent strands)
  (u.2 = v.2 ∧ (u.1.val + 1 = v.1.val ∨ v.1.val + 1 = u.1.val)) ∨
  -- Crossing edges
  (∃ (t : Fin b.length),
     let i := (b.get t).val
     ((u.1.val = i ∧ u.2 = t.castSucc ∧ v.1.val = i + 1 ∧ v.2 = t.succ) ∨
      (u.1.val = i + 1 ∧ u.2 = t.castSucc ∧ v.1.val = i ∧ v.2 = t.succ)))

instance {n : ℕ} (b : Braid n) : DecidableRel (AugRhomboidAdj b) := by
  intro u v; unfold AugRhomboidAdj; exact inferInstance

/-- The augmented rhomboid graph of a braid. -/
def AugmentedRhomboidGraph {n : ℕ} (b : Braid n) : SimpleGraph (Fin n × Fin (b.length + 1)) :=
  SimpleGraph.fromRel (AugRhomboidAdj b)

instance augDecAdj {n : ℕ} (b : Braid n) : DecidableRel (AugmentedRhomboidGraph b).Adj := by
  intro u v; simp only [AugmentedRhomboidGraph, SimpleGraph.fromRel_adj]; exact inferInstance

/-! ## Abstract Simplicial Complex -/

/-- An abstract simplicial complex on vertex type V. -/
structure AbstractSimplicialComplex (V : Type*) where
  faces : Set (Finset V)
  down_closed : ∀ f ∈ faces, ∀ g ⊆ f, g ∈ faces

/-- The independence complex of a simple graph. -/
def IndependenceComplex {V : Type*} [DecidableEq V] (G : SimpleGraph V) :
    AbstractSimplicialComplex V where
  faces := {s : Finset V | ∀ u ∈ s, ∀ v ∈ s, u ≠ v → ¬G.Adj u v}
  down_closed := by
    intro f hf g hg u hu v hv huv
    exact hf u (hg hu) v (hg hv) huv

/-! ## Free Pairs and Collapse -/

variable {V : Type*} [DecidableEq V]

/-- A free pair (σ, τ) in K: σ ⊂ τ with |τ| = |σ| + 1, and σ is a face of no
    other face besides σ and τ. -/
def is_free_pair (K : AbstractSimplicialComplex V) (σ τ : Finset V) : Prop :=
  σ ∈ K.faces ∧ τ ∈ K.faces ∧
  σ ⊆ τ ∧ τ.card = σ.card + 1 ∧
  (∀ ρ ∈ K.faces, σ ⊆ ρ → ρ = σ ∨ ρ = τ)

/-
Elementary collapse: remove σ and τ from the complex.
-/
def elementary_collapse (K : AbstractSimplicialComplex V) (σ τ : Finset V)
    (h : is_free_pair K σ τ) : AbstractSimplicialComplex V where
  faces := K.faces \ {σ, τ}
  down_closed := by
    simp_all +decide [ is_free_pair ];
    intro f hf hσ hτ g hg
    have hg_in_K : g ∈ K.faces := by
      exact K.down_closed f hf g hg
    have hg_ne_σ : g ≠ σ := by
      grind +ring
    have hg_ne_τ : g ≠ τ := by
      grind +ring
    exact ⟨hg_in_K, hg_ne_σ, hg_ne_τ⟩

/-- The collapse relation: K collapses to L via a sequence of elementary collapses. -/
inductive Collapses : AbstractSimplicialComplex V → AbstractSimplicialComplex V → Prop where
  | refl (K : AbstractSimplicialComplex V) : Collapses K K
  | step {K L : AbstractSimplicialComplex V} (σ τ : Finset V) (h : is_free_pair L σ τ) :
      Collapses K L → Collapses K (elementary_collapse L σ τ h)

lemma Collapses_trans {K L M : AbstractSimplicialComplex V} :
    Collapses K L → Collapses L M → Collapses K M := by
  intro hKL hLM
  induction hLM with
  | refl => exact hKL
  | step σ τ h _ ih => exact Collapses.step σ τ h ih

/-! ## Connectivity -/

/-- Connectivity of a simplicial complex: any two vertices (singletons in the complex)
    are connected by a path of edges (2-element faces). -/
def is_connected (K : AbstractSimplicialComplex V) : Prop :=
  ∀ u v, {u} ∈ K.faces → {v} ∈ K.faces →
  ∃ (p : List V), p ≠ [] ∧ p.head? = some u ∧ p.getLast? = some v ∧
    List.Chain' (fun a b => {a, b} ∈ K.faces) p

/-- A complex has no free pairs. -/
def has_no_free_pairs (K : AbstractSimplicialComplex V) : Prop :=
  ∀ σ τ, ¬is_free_pair K σ τ
