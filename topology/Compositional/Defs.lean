/-
  Compositional Fixed Strand Conjecture: Definitions

  This file defines:
  - Braid composition (horizontal composition = list append)
  - Interface/boundary vertex sets
  - Regional classification of simplices and free pairs
  - Boundary states of collapsed complexes
  - The abstract Khovanov framework
-/
import Mathlib
import FixedStrand.Defs
import FixedStrand.Collapse

set_option maxHeartbeats 800000

open scoped Classical

/-! ## Braid Composition -/

/-- Horizontal composition of two braids on n strands is list concatenation. -/
def braid_compose {n : ℕ} (b1 b2 : Braid n) : Braid n := by
  unfold Braid at *; exact b1 ++ b2

lemma braid_compose_length {n : ℕ} (b1 b2 : Braid n) :
    (braid_compose b1 b2).length = b1.length + b2.length := by
  simp [braid_compose, Braid]

/-! ## Vertex Regions in the Composed Graph -/

/-- The vertex type of the augmented rhomboid graph of a composed braid. -/
abbrev ComposedVertex (n : ℕ) (b1 b2 : Braid n) :=
  Fin n × Fin ((braid_compose b1 b2).length + 1)

/-- The interface time step where the two braids are glued. -/
def interface_time {n : ℕ} (b1 : Braid n) : ℕ := b1.length

/-- A vertex is in the B1 region (left part) if its time coordinate is ≤ the interface. -/
def in_B1_region {n : ℕ} (b1 b2 : Braid n) (v : ComposedVertex n b1 b2) : Prop :=
  v.2.val ≤ interface_time b1

/-- A vertex is in the B2 region (right part) if its time coordinate is ≥ the interface. -/
def in_B2_region {n : ℕ} (b1 b2 : Braid n) (v : ComposedVertex n b1 b2) : Prop :=
  v.2.val ≥ interface_time b1

/-- A vertex is at the interface if its time coordinate equals the interface time. -/
def at_interface {n : ℕ} (b1 b2 : Braid n) (v : ComposedVertex n b1 b2) : Prop :=
  v.2.val = interface_time b1

/-- A vertex is in the strict interior of B1 if its time is strictly less than the interface. -/
def in_strict_B1 {n : ℕ} (b1 b2 : Braid n) (v : ComposedVertex n b1 b2) : Prop :=
  v.2.val < interface_time b1

/-- A vertex is in the strict interior of B2 if its time is strictly greater than the interface. -/
def in_strict_B2 {n : ℕ} (b1 b2 : Braid n) (v : ComposedVertex n b1 b2) : Prop :=
  v.2.val > interface_time b1

/-- A vertex is in the interface neighborhood (within 1 time step of the interface). -/
def in_interface_nbhd {n : ℕ} (b1 b2 : Braid n) (v : ComposedVertex n b1 b2) : Prop :=
  v.2.val + 1 ≥ interface_time b1 ∧ interface_time b1 + 1 ≥ v.2.val

/-! ## Regional Classification of Simplices -/

/-- A simplex (finite set of vertices) is B1-interior if all its vertices are in strict B1. -/
def simplex_B1_interior {n : ℕ} (b1 b2 : Braid n)
    (s : Finset (ComposedVertex n b1 b2)) : Prop :=
  ∀ v ∈ s, in_strict_B1 b1 b2 v

/-- A simplex is B2-interior if all its vertices are in strict B2. -/
def simplex_B2_interior {n : ℕ} (b1 b2 : Braid n)
    (s : Finset (ComposedVertex n b1 b2)) : Prop :=
  ∀ v ∈ s, in_strict_B2 b1 b2 v

/-- A simplex touches the interface if at least one vertex is in the interface neighborhood. -/
def simplex_touches_interface {n : ℕ} (b1 b2 : Braid n)
    (s : Finset (ComposedVertex n b1 b2)) : Prop :=
  ∃ v ∈ s, in_interface_nbhd b1 b2 v

/-! ## Interface Free Pairs -/

/-- A free pair is interface-local if σ (the smaller face) contains a vertex
    in the interface neighborhood. This ensures that any face ρ ⊇ σ also
    touches the interface (since σ ⊆ ρ), which is crucial for the locality proof. -/
def is_interface_free_pair {n : ℕ} (b1 b2 : Braid n)
    (K : AbstractSimplicialComplex (ComposedVertex n b1 b2))
    (σ τ : Finset (ComposedVertex n b1 b2)) : Prop :=
  is_free_pair K σ τ ∧ (∃ v ∈ σ, in_interface_nbhd b1 b2 v)

/-! ## Interface Subgraph -/

/-- The interface subgraph consists of vertices within 1 time step of the gluing point. -/
def InterfaceSubgraph {n : ℕ} (b1 b2 : Braid n) :
    SimpleGraph { v : ComposedVertex n b1 b2 // in_interface_nbhd b1 b2 v } :=
  (AugmentedRhomboidGraph (braid_compose b1 b2)).comap Subtype.val

/-! ## Boundary State -/

/-- The boundary state of a simplicial complex at a time step t is the restriction of the
    complex to vertices at time t, i.e., the set of faces consisting only of vertices
    at time t. This captures the "outgoing data" of a collapsed complex at the boundary. -/
def boundary_faces_at {n : ℕ} {b : Braid n}
    (K : AbstractSimplicialComplex (Fin n × Fin (b.length + 1)))
    (t : ℕ) : Set (Finset (Fin n × Fin (b.length + 1))) :=
  { s ∈ K.faces | ∀ v ∈ s, v.2.val = t }

/-- Two complexes have the same boundary state at time t if they have the same faces at t. -/
def same_boundary_at {n : ℕ} {b : Braid n}
    (K L : AbstractSimplicialComplex (Fin n × Fin (b.length + 1)))
    (t : ℕ) : Prop :=
  boundary_faces_at K t = boundary_faces_at L t

/-! ## Abstract Khovanov Framework -/

/-- An abstract Khovanov framework packages the algebraic structures needed for the
    Compositional FSC. This captures the key properties of Khovanov homology without
    requiring a full formalization of TQFTs.

    The original FSC (already proved) establishes that free pairs in the independence
    complex correspond to algebraic reductions. This framework axiomatizes that
    correspondence so we can extend it to compositions. -/
class KhovanovFramework (n : ℕ) where
  /-- The type of (reduced) Khovanov chain complexes. -/
  KhComplex : Type
  /-- The minimal Khovanov complex of a braid. -/
  khovanov_min : Braid n → KhComplex
  /-- The unreduced Khovanov bracket (before Morse collapse). -/
  khovanov_full : Braid n → KhComplex
  /-- Tensor product of two complexes. -/
  tensor : KhComplex → KhComplex → KhComplex
  /-- The type of algebraic isomorphisms (saddle isomorphisms / Gaussian eliminations). -/
  AlgReduction : Type
  /-- The set of algebraic reductions needed to go from unreduced to reduced. -/
  full_reductions : Braid n → Set AlgReduction
  /-- The algebraic reductions from a tensor product to the minimal complex. -/
  tensor_reductions : KhComplex → KhComplex → Set AlgReduction
  /-- The FSC correspondence: free pairs ↔ algebraic reductions (the original FSC). -/
  fsc_bij : ∀ b : Braid n,
    ∃ f : { p : Finset (Fin n × Fin (b.length + 1)) × Finset (Fin n × Fin (b.length + 1)) //
            is_free_pair (IndependenceComplex (AugmentedRhomboidGraph b)) p.1 p.2 } →
           full_reductions b,
      Function.Bijective f
  /-- The reductions for a composed braid decompose into parts + interface. -/
  reduction_decomposition : ∀ b1 b2 : Braid n,
    ∃ (R1 : Set AlgReduction) (R2 : Set AlgReduction) (R_glue : Set AlgReduction),
      full_reductions (braid_compose b1 b2) = R1 ∪ R2 ∪ R_glue ∧
      R_glue = tensor_reductions (khovanov_min b1) (khovanov_min b2)
  /-- The FSC bijection maps interface free pairs to tensor (interface) reductions.
      This captures the key structural property: the discrete Morse matching
      constructed in the FSC proof is compatible with the regional decomposition
      of the composed graph. Interface free pairs correspond to the "new"
      algebraic reductions needed when forming the tensor product. -/
  fsc_interface_bij : ∀ b1 b2 : Braid n,
    ∃ (f : { p : Finset (ComposedVertex n b1 b2) × Finset (ComposedVertex n b1 b2) //
             is_interface_free_pair b1 b2
               (IndependenceComplex (AugmentedRhomboidGraph (braid_compose b1 b2))) p.1 p.2 } →
           tensor_reductions (khovanov_min b1) (khovanov_min b2)),
      Function.Bijective f

/-! ## The Glue Matching -/

/-- The glue matching M_glue: the set of interface-local free pairs in the composed
    independence complex. These are the free pairs whose smaller face σ contains
    a vertex in the interface neighborhood. -/
def glue_matching {n : ℕ} (b1 b2 : Braid n) :
    Set (Finset (ComposedVertex n b1 b2) × Finset (ComposedVertex n b1 b2)) :=
  { p | is_interface_free_pair b1 b2
        (IndependenceComplex (AugmentedRhomboidGraph (braid_compose b1 b2))) p.1 p.2 }
