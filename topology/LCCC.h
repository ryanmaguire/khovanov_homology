/*
 *  LCCC.h
 *
 *  Linear Combination of Canned Cobordisms with coefficients in either
 *  Z or a prime field F_p.
 *
 *  Purpose:
 *      Provides the concrete algebraic coefficient type used in the
 *      Khovanov chain complex. Coefficients remain stored as int values,
 *      but arithmetic is interpreted according to one active coefficient
 *      mode selected before the complex is constructed.
 *
 *  Structure:
 *      An LCCC is a singly-linked list of (coefficient, cobordism)
 *      pairs. Structurally equal cobordisms are combined by adding
 *      coefficients, and zero coefficients are removed.
 *
 *  Coefficient modes:
 *      KH_COEFF_Z   : ordinary integer arithmetic. Units are +/-1.
 *      KH_COEFF_FP  : arithmetic modulo a prime p. Every nonzero
 *                     coefficient is a unit.
 *
 *  Important:
 *      Select the coefficient mode once, before constructing a Komplex,
 *      and do not change it while LCCC objects from that computation are
 *      still in use.
 */

#ifndef LCCC_H
#define LCCC_H

#include "CannedCobordism.h"
#include <stdbool.h>

/* ================================================================
 *  Coefficient mode
 * ================================================================ */

typedef enum KHCoefficientMode {
  KH_COEFF_Z = 0,
  KH_COEFF_FP = 1
} KHCoefficientMode;

/*
 * Restore ordinary integer coefficients. This is also the default mode.
 */
void LCCC_setCoefficientIntegers(void);

/*
 * Use coefficients in F_p. Returns false unless p is prime.
 * Call this before constructing the Khovanov complex.
 */
bool LCCC_setCoefficientModPrime(int p);

KHCoefficientMode LCCC_getCoefficientMode(void);
int LCCC_getCoefficientModulus(void);

/*
 * Small coefficient helpers exposed so Komplex/TangleKomplex can make the
 * same unit/sign decisions as LCCC without duplicating arithmetic rules.
 */
int LCCC_coeffNormalize(int value);
int LCCC_coeffAdd(int a, int b);
int LCCC_coeffMultiply(int a, int b);
int LCCC_coeffNegate(int value);
bool LCCC_coeffIsZero(int value);
bool LCCC_coeffIsUnit(int value);
bool LCCC_coeffInverse(int value, int *inverse_out);

/* ================================================================
 *  Type Definitions
 * ================================================================ */

/*
 *  LCCCTerm — a single node in the linked list.
 *
 *  coeff is always an int. In KH_COEFF_Z it is an ordinary integer;
 *  in KH_COEFF_FP it is stored in canonical form 0,...,p-1.
 */
typedef struct LCCCTerm {
  int coeff;
  CannedCobordism *cobordism;
  struct LCCCTerm *next;
} LCCCTerm;

/*
 *  LCCC — the linear combination itself.
 */
typedef struct LCCC {
  LCCCTerm *head;
  int count; /* Number of nonzero terms, for quick zero checks. */
} LCCC;

/*
 *  RingElement — lightweight scalar wrapper retained for API compatibility
 *  with CobMatrix. Its value is interpreted using the active coefficient
 *  mode above.
 */
typedef struct RingElement {
  int value;
} RingElement;

/* ================================================================
 *  Construction / Destruction
 * ================================================================ */

/* Create an empty LCCC (the zero element). */
LCCC *LCCC_createZero(void);

/*
 * Create an LCCC containing a single cobordism term. The coefficient is
 * normalized according to the active coefficient mode; a zero coefficient
 * produces the zero LCCC.
 */
LCCC *LCCC_createSingle(CannedCobordism *cc, int coeff);

/* Deep-clone an LCCC list; cobordism pointers remain shared. */
LCCC *LCCC_clone(const LCCC *lc);

/* Free an LCCC and its term list, but not the cobordism pointers. */
void LCCC_free(LCCC *lc);

/* ================================================================
 *  Arithmetic
 * ================================================================ */

/* Add two LCCCs using the active coefficient arithmetic. */
LCCC *LCCC_add(LCCC *a, LCCC *b);

/*
 * Compose two LCCCs, distributing composition over all term pairs and
 * multiplying their coefficients in the active coefficient system.
 */
LCCC *LCCC_compose(LCCC *a, LCCC *b);

/* Multiply every term by a scalar interpreted in the active mode. */
LCCC *LCCC_multiply(LCCC *a, RingElement *coeff);

/*
 * Reduce using the Bar-Natan surface relations. Numerical factors such as
 * torus = 2 are interpreted in the active coefficient system, so for example
 * torus = 0 over F_2.
 */
LCCC *LCCC_reduce(LCCC *a);

/*
 * Invert a single-term unit pivot. Over Z only +/-1 are units. Over F_p
 * every nonzero coefficient is a unit and is replaced by its modular inverse.
 * The current Gaussian-elimination path uses identity-like cobordism pivots,
 * whose topological inverse is represented by the same cobordism term.
 */
LCCC *LCCC_invert(LCCC *lc);

/* Negate an LCCC in the active coefficient system. */
LCCC *LCCC_negate(LCCC *lc);

LCCC *LCCC_capOffTop(LCCC *lc, Cap *new_top, bool add_dot);
LCCC *LCCC_cupOnBottom(LCCC *lc, Cap *new_bottom, bool add_dot);

/* ================================================================
 *  Predicates
 * ================================================================ */

bool LCCC_isZero(LCCC *lc);
bool LCCC_contains(const LCCC *lc, const CannedCobordism *cc);

#endif // LCCC_H
