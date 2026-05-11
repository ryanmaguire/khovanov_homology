#ifndef COBMATRIX_H
#define COBMATRIX_H

#include <stdbool.h>
#include <stdlib.h>

// Forward declarations for topological / algebraic types
// In the C version, LCCC (Linear Combination of Canned Cobordisms) and Caps
// are represented abstractly.
typedef struct LCCC LCCC; 
typedef struct Cap Cap;
typedef struct RingElement RingElement;

LCCC* LCCC_add(LCCC* a, LCCC* b);
LCCC* LCCC_compose(LCCC* a, LCCC* b);
LCCC* LCCC_multiply(LCCC* a, RingElement* coeff);
LCCC* LCCC_reduce(LCCC* a);
bool LCCC_isZero(LCCC* a);
void LCCC_free(LCCC* a);

// Core Data Structures
struct MatrixEntry {
    int column_index;
    LCCC* value;
    struct MatrixEntry* next;
};

typedef struct MatrixEntry MatrixEntry;

typedef struct MatrixRow {
    MatrixEntry* head;
} MatrixRow;

typedef struct SmoothingColumn {
    int n;
    int* numbers; // Array of integers
    Cap** smoothings; // Array of Cap pointers
} SmoothingColumn;

SmoothingColumn* SmoothingColumn_clone(SmoothingColumn* col);
bool SmoothingColumn_equals(SmoothingColumn* a, SmoothingColumn* b);

typedef struct CobMatrix {
    SmoothingColumn* source;
    SmoothingColumn* target;
    MatrixRow* entries; // Array of MatrixRow of size 'target->n'
} CobMatrix;

// Constructor & Destructor
CobMatrix* CobMatrix_create(SmoothingColumn* source, SmoothingColumn* target, bool shared);
void CobMatrix_free(CobMatrix* m);

// Core Matrix Operations
void CobMatrix_putEntry(CobMatrix* m, int row_idx, int col_idx, LCCC* lc);
void CobMatrix_addEntry(CobMatrix* m, int row_idx, int col_idx, LCCC* lc);
LCCC** CobMatrix_unpackRow(CobMatrix* m, int row_idx);

// Mathematical Operations
CobMatrix* CobMatrix_compose(CobMatrix* this_m, CobMatrix* that_m);
void CobMatrix_multiply(CobMatrix* m, RingElement* n);
void CobMatrix_add(CobMatrix* dest, CobMatrix* src);
void CobMatrix_reduce(CobMatrix* m);

// Properties / Checks
bool CobMatrix_isZero(CobMatrix* m);
bool CobMatrix_equals(CobMatrix* a, CobMatrix* b);

// Row / Column Extraction
// Modifies 'm' by removing the row/column, and returns a new CobMatrix representing the slice.
CobMatrix* CobMatrix_extractColumn(CobMatrix* m, int column_idx);
CobMatrix* CobMatrix_extractRow(CobMatrix* m, int row_idx);

#endif // COBMATRIX_H
