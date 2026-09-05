#ifndef INT_MATRIX_H
#define INT_MATRIX_H

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct Mat Mat;

struct Mat {
    int rows;
    int cols;
    int64_t** matrix; /* replacement of Java's BigInteger[][] */

    /* chaining complexes */
    Mat* prev;
    Mat* next;

    /* tracking arrays — replacement for Java's List<Integer> source, target */
    int* source;
    int* target;
};

/* memory management */
Mat* createMat(int rows, int cols);
void freeMat(Mat* m);
void isolateMat(Mat* m);

/* status checks */
bool isDiag(Mat* m);
bool isnull(Mat* m);

/* gaussian elim and smith normal — column ops */
void multColumn(Mat* m, int col, int64_t scalar);
void SwapCols(Mat* m, int col1, int col2);
void addColumn(Mat* m, int dest_col, int src_col, int64_t scalar);

/* row operations */
void multRow(Mat* m, int row, int64_t scalar);
void SwapRows(Mat* m, int row1, int row2);
void addRow(Mat* m, int dest_row, int src_row, int64_t scalar);

/* other functions */
int rowNonZeroes(Mat* m, int row);
int columnNonZeroes(Mat* m, int col);
int zeroRowsToEnd(Mat* m);
int zeroColumnsToEnd(Mat* m);

void toSmithForm(Mat* m);
void printMat(Mat* m);

/*
 * Extract free image-rank and torsion coefficients from a matrix that is
 * already in Smith normal form (call toSmithForm first).
 *
 * Semantics (matches the enhanced-state / Khovanov pipeline):
 *   - free_rank  = number of non-zero diagonal entries
 *                  = rank of the image of this differential over Q
 *   - torsions[k] = |d_ii| whenever |d_ii| > 1
 *                  (invariant factors that produce Z/d torsion)
 *
 * Parameters:
 *   m          – matrix in SNF (not modified)
 *   free_rank  – out; may be NULL if not needed
 *   torsions   – out buffer of length at least max_tors; may be NULL
 *   max_tors   – capacity of torsions[]
 *   n_tors     – out; number of torsion coefficients written (capped at max_tors)
 *
 * Returns 0 on success, -1 if m is NULL.
 *
 * Note: Betti numbers of the chain complex are recovered upstream as
 *   β_h = dim C_h − rank(d_h) − rank(d_{h-1})
 * This helper only reports rank(d) and the torsion part of coker(d).
 */
int smith_extract(const Mat* m,
                  int* free_rank,
                  int* torsions,
                  int max_tors,
                  int* n_tors);

/*
 * Convenience: run toSmithForm then smith_extract.
 * Mutates m into SNF. Same out-parameter contract as smith_extract.
 */
int toSmithForm_extract(Mat* m,
                        int* free_rank,
                        int* torsions,
                        int max_tors,
                        int* n_tors);

#endif /* INT_MATRIX_H */
