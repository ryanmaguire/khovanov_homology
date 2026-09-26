#ifndef DTTOPD_H
#define DTTOPD_H

#include <stdbool.h>
#include <stddef.h>

/* Forward declaration from PDScanner.h. */
typedef struct PDDiagram PDDiagram;

/*
 * This converter is intended for the compact alphabetic DT data used by the
 * knot-table pipeline (<= 19 crossings).  The alphabet itself can encode up
 * to 26 crossings, but the current planar-embedding search is deliberately
 * capped at 20 so its worst-case search remains small and predictable.
 */
#define DTTOPD_MAX_CROSSINGS 20

/*
 * Decode compact alphabetical Dowker-Thistlethwaite notation.
 *
 *   a,b,c,...  ->  +2,+4,+6,...
 *   A,B,C,...  ->  -2,-4,-6,...
 *
 * ASCII whitespace is ignored.  On success, *dt_out is heap allocated and
 * must be freed by the caller with free().
 */
bool DTtoPD_decodeAlphabetical(const char *code,
                               int **dt_out,
                               int *crossing_count_out,
                               char *reason,
                               size_t reason_size);

/*
 * Convert a numerical DT code directly into the PDDiagram representation used
 * by PDScanner.  The resulting PD edge labels are already dense and 0-based:
 * 0,1,...,2n-1.  diagram->signs is filled at the same time.
 *
 * DT notation has an unavoidable global mirror ambiguity.  This routine makes
 * the choice deterministic by selecting the planar rotation system whose first
 * crossing has sign +1.  The globally mirrored system has all signs reversed.
 */
bool DTtoPD_fromNumeric(PDDiagram *diagram,
                        const int *dt,
                        int crossing_count,
                        char *reason,
                        size_t reason_size);

/*
 * Convenience wrapper: decode an alphabetical DT code and convert it directly
 * into a PDDiagram suitable for PDScanner_build.  Equivalent to
 * DTtoPD_decodeAlphabetical followed by DTtoPD_fromNumeric.
 */
bool DTtoPD_fromAlphabetical(PDDiagram *diagram,
                             const char *code,
                             char *reason,
                             size_t reason_size);

#endif /* DTTOPD_H */
