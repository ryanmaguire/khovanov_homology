/*
 * Purpose:
 *      The identity cobordism: source and target are the same Cap,
 *      composition returns the other cobordism unchanged, and
 *      isIsomorphism is always true.
 */

#ifndef IDENTITYCANNEDCOBORDISM_H
#define IDENTITYCANNEDCOBORDISM_H

#include "CannedCobordism.h"

/*
 * Purpose:
 *      Creates an identity CannedCobordism for the given Cap.
 *      The returned struct has its vtable function pointers populated
 *      with identity-specific implementations.
 * Arguments:
 *      cap
 * Argument descriptions:
 *      cap - The Cap to use as both source and target.
 * Output:
 *      CannedCobordism pointer
 * Output description:
 *      Returns a newly allocated CannedCobordism. The caller takes ownership
 *      of the returned pointer. Free with CannedCobordism_free().
 * Method:
 *      Allocates a CannedCobordism, sets its source and target to the given
 *      Cap, and wires its vtable pointers to the local identity functions.
 */
CannedCobordism *IdentityCannedCobordism_create(Cap *cap);

#endif /* IDENTITYCANNEDCOBORDISM_H */
