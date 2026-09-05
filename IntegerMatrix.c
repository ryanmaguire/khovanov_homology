#include "IntegerMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static void integerMatrixFatal(const char* message) {
    fprintf(stderr, "FATAL: IntegerMatrix: %s\n", message);
    abort();
}

static int64_t checkedAdd(int64_t a, int64_t b) {
    if ((b > 0 && a > INT64_MAX - b) ||
        (b < 0 && a < INT64_MIN - b)) {
        integerMatrixFatal("int64_t addition overflow");
    }
    return a + b;
}

static int64_t checkedMul(int64_t a, int64_t b) {
    if (a == 0 || b == 0) return 0;

    if (a > 0) {
        if (b > 0) {
            if (a > INT64_MAX / b) {
                integerMatrixFatal("int64_t multiplication overflow");
            }
        } else {
            if (b < INT64_MIN / a) {
                integerMatrixFatal("int64_t multiplication overflow");
            }
        }
    } else {
        if (b > 0) {
            if (a < INT64_MIN / b) {
                integerMatrixFatal("int64_t multiplication overflow");
            }
        } else {
            if (b < INT64_MAX / a) {
                integerMatrixFatal("int64_t multiplication overflow");
            }
        }
    }

    return a * b;
}

static int64_t i64_abs(int64_t v) {
    return (v < 0) ? -v : v;
}

/* * Allocates memory for the matrix structure and its internal 2D grid.
 * Initializes tracking pointers and links to NULL.
 */
Mat* createMat(int rows, int cols) {
    if (rows < 0 || cols < 0) return NULL;

    Mat* m = (Mat*)malloc(sizeof(Mat));
    if (!m) return NULL;

    m->rows = rows;
    m->cols = cols;

    /* Allocate array of pointers for rows using 64-bit type */
    m->matrix = (int64_t**)malloc((size_t)rows * sizeof(int64_t*));
    if (rows > 0 && !m->matrix) {
        free(m);
        return NULL;
    }

    for (int i = 0; i < rows; i++) {
        /* Allocate and zero-initialize columns for each row */
        m->matrix[i] = (int64_t*)calloc((size_t)cols, sizeof(int64_t));
        if (cols > 0 && !m->matrix[i]) {
            for (int j = 0; j < i; j++) {
                free(m->matrix[j]);
            }
            free(m->matrix);
            free(m);
            return NULL;
        }
    }

    m->prev = NULL;
    m->next = NULL;
    m->source = NULL;
    m->target = NULL;
    return m;
}

/* * Safely deallocates the matrix and its associated tracking arrays.
 */
void freeMat(Mat* m) {
    if (!m) return;

    /* Free each individual row row-by-row */
    for (int i = 0; i < m->rows; i++) {
        free(m->matrix[i]);
    }
    free(m->matrix);

    /* Free generator tracking arrays if they were allocated */
    if (m->source != NULL) free(m->source);
    if (m->target != NULL) free(m->target);

    free(m);
}

/* * Isolates the matrix by breaking its links within the chain complex.
 */
void isolateMat(Mat* m) {
    if (!m) return;

    /* Safely disconnect from the neighbors in the chain complex */
    if (m->prev != NULL) m->prev->next = NULL;
    if (m->next != NULL) m->next->prev = NULL;

    m->prev = NULL;
    m->next = NULL;
}

/* Checks if the matrix is strictly diagonal */
bool isDiag(Mat* m) {
    int rows = m->rows, cols = m->cols;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (i != j && m->matrix[i][j] != 0) {
                return false;
            }
        }
    }
    return true;
}

/* Checks if the matrix is completely zero (null matrix) */
bool isnull(Mat* m) {
    int rows = m->rows, cols = m->cols;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (m->matrix[i][j] != 0) {
                return false;
            }
        }
    }
    return true;
}

/* Counts non-zero elements in a specific row */
int rowNonZeroes(Mat* m, int row) {
    int count = 0;
    for (int j = 0; j < m->cols; j++) {
        if (m->matrix[row][j] != 0) count++;
    }
    return count;
}

/* Counts non-zero elements in a specific column */
int columnNonZeroes(Mat* m, int col) {
    int count = 0;
    for (int i = 0; i < m->rows; i++) {
        if (m->matrix[i][col] != 0) count++;
    }
    return count;
}

/* Low-level operation: Swaps two row pointers in the current matrix (O(1) complexity) */
void swaprows(Mat* m, int row1, int row2) {
    int64_t* temp = m->matrix[row1];
    m->matrix[row1] = m->matrix[row2];
    m->matrix[row2] = temp;
}

/* Low-level operation: Swaps two columns element by element (O(Rows) complexity) */
void swapcols(Mat* m, int col1, int col2) {
    for (int i = 0; i < m->rows; i++) {
        int64_t temp = m->matrix[i][col1];
        m->matrix[i][col1] = m->matrix[i][col2];
        m->matrix[i][col2] = temp;
    }
}

/* Low-level operation: Row Addition -> Row(a) = Row(a) + n * Row(b) */
void addRow2(Mat* m, int a, int b, int64_t n) {
    for (int i = 0; i < m->cols; i++) {
        int64_t product = checkedMul(m->matrix[b][i], n);
        m->matrix[a][i] = checkedAdd(m->matrix[a][i], product);
    }
}

/* Low-level operation: Column Addition -> Col(a) = Col(a) + n * Col(b) */
void addColumn2(Mat* m, int a, int b, int64_t n) {
    for (int i = 0; i < m->rows; i++) {
        int64_t product = checkedMul(m->matrix[i][b], n);
        m->matrix[i][a] = checkedAdd(m->matrix[i][a], product);
    }
}

/* Low-level operation: Multiplies a row by a scalar factor */
void multRow2(Mat* m, int a, int64_t n) {
    if (n != 1 && n != -1) {
        integerMatrixFatal("non-unimodular row scaling");
    }

    for (int i = 0; i < m->cols; i++) {
        m->matrix[a][i] = checkedMul(m->matrix[a][i], n);
    }
}

/* Low-level operation: Multiplies a column by a scalar factor */
void multColumn2(Mat* m, int a, int64_t n) {
    if (n != 1 && n != -1) {
        integerMatrixFatal("non-unimodular column scaling");
    }

    for (int i = 0; i < m->rows; i++) {
        m->matrix[i][a] = checkedMul(m->matrix[i][a], n);
    }
}

/* * High-level operation: Cascading Row Swap.
 * Changing the basis of the target space of d_n (swapping rows) requires
 * swapping the columns of d_{n+1} to maintain the chain complex property (d ◦ d = 0).
 */
void SwapRows(Mat* m, int row1, int row2) {
    if (row1 == row2) return;
    swaprows(m, row1, row2);

    /* Cascade to the next operator's columns to maintain algebraic consistency */
    if (m->next != NULL) {
        swapcols(m->next, row1, row2);
    }

    /* Update target generator indices mappings */
    if (m->target != NULL) {
        int temp = m->target[row1];
        m->target[row1] = m->target[row2];
        m->target[row2] = temp;
    }
}

/* * High-level operation: Cascading Column Swap.
 * Swapping columns in d_n changes the domain basis, requiring a corresponding
 * row swap in d_{n-1} to maintain the chain complex equation.
 */
void SwapCols(Mat* m, int col1, int col2) {
    if (col1 == col2) return;
    swapcols(m, col1, col2);

    /* Cascade to the previous operator's rows */
    if (m->prev != NULL) {
        swaprows(m->prev, col1, col2);
    }

    /* Update source generator indices mappings */
    if (m->source != NULL) {
        int temp = m->source[col1];
        m->source[col1] = m->source[col2];
        m->source[col2] = temp;
    }
}

/* Interfaces and wrappers for elementary operations */
void addRow(Mat* m, int a, int b, int64_t n) { addRow2(m, a, b, n); }
void addColumn(Mat* m, int a, int b, int64_t n) { addColumn2(m, a, b, n); }
void multRow(Mat* m, int a, int64_t n) { multRow2(m, a, n); }
void multColumn(Mat* m, int a, int64_t n) { multColumn2(m, a, n); }

/* Shifts all rows that contain only zeroes to the bottom of the matrix */
int zeroRowsToEnd(Mat* m) {
    int nzrows = m->rows;
    for (int i = 0; i < nzrows; i++) {
        while (rowNonZeroes(m, i) == 0 && i < nzrows) {
            SwapRows(m, i, --nzrows);
        }
    }
    return nzrows;
}

/* Shifts all columns that contain only zeroes to the rightmost side of the matrix */
int zeroColumnsToEnd(Mat* m) {
    int nzcols = m->cols;
    for (int i = 0; i < nzcols; i++) {
        while (columnNonZeroes(m, i) == 0 && i < nzcols) {
            SwapCols(m, i, --nzcols);
        }
    }
    return nzcols;
}

/* * Computes the Smith Normal Form (SNF) of the matrix over the Ring of Integers (Z).
 * This executes Euclidean reduction to reveal the free ranks (Betti numbers)
 * and structural torsion coefficients of the Khovanov Homology group.
 */
void toSmithForm(Mat* m) {
    isolateMat(m);

    int row = 0, col = 0;
    while (row < m->rows && col < m->cols) {
        while (row < m->rows && rowNonZeroes(m, row) == 0) row++;
        while (col < m->cols && columnNonZeroes(m, col) == 0) col++;

        if (row >= m->rows || col >= m->cols) break;

        if (row > col) { SwapRows(m, row, col); row = col; }
        else if (col > row) { SwapCols(m, row, col); col = row; }

        /* Euclidean elimination loop */
        while (rowNonZeroes(m, row) > 1 || columnNonZeroes(m, col) > 1 || m->matrix[row][col] <= 0) {

            if (rowNonZeroes(m, row) == 0 && columnNonZeroes(m, col) == 0) break;

            /* clear column */
            while (columnNonZeroes(m, col) > 1 || m->matrix[row][col] <= 0) {
                /* force positives */
                for (int j = row; j < m->rows; j++) {
                    if (m->matrix[j][col] < 0) multRow(m, j, -1);
                }

                int64_t min_val = -1;
                int idxmin = -1;
                for (int j = row; j < m->rows; j++) {
                    int64_t v = m->matrix[j][col];
                    if (v > 0 && (min_val == -1 || v < min_val)) {
                        min_val = v; idxmin = j;
                    }
                }

                if (min_val == -1) break;
                if (idxmin != row) SwapRows(m, row, idxmin);

                for (int j = row + 1; j < m->rows; j++) {
                    if (m->matrix[j][col] != 0) {
                        addRow(m, j, row, -(m->matrix[j][col] / min_val));
                    }
                }
            }

            /* clear the row */
            while (rowNonZeroes(m, row) > 1 || m->matrix[row][col] <= 0) {
                /* force positive values */
                for (int j = col; j < m->cols; j++) {
                    if (m->matrix[row][j] < 0) multColumn(m, j, -1);
                }

                int64_t min_val = -1;
                int idxmin = -1;
                for (int j = col; j < m->cols; j++) {
                    int64_t v = m->matrix[row][j];
                    if (v > 0 && (min_val == -1 || v < min_val)) {
                        min_val = v; idxmin = j;
                    }
                }

                if (min_val == -1) break;
                if (idxmin != col) SwapCols(m, col, idxmin);

                for (int j = col + 1; j < m->cols; j++) {
                    if (m->matrix[row][j] != 0) {
                        addColumn(m, j, col, -(m->matrix[row][j] / min_val));
                    }
                }
            }
        }

        /* move along diag */
        row++;
        col++;
    }

    zeroRowsToEnd(m);
    zeroColumnsToEnd(m);
}

/* Outputs the matrix values onto stdout with standard formatting */
void printMat(Mat* m) {
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            printf("%5lld", (long long)m->matrix[i][j]);
        }
        printf("\n");
    }
}

/*
 * Read free image-rank and torsion coefficients from an already-SNF matrix.
 * See IntegerMatrix.h for full contract.
 */
int smith_extract(const Mat* m,
                  int* free_rank,
                  int* torsions,
                  int max_tors,
                  int* n_tors)
{
    if (!m) return -1;

    int rank = 0;
    int nt = 0;
    int limit = (m->rows < m->cols) ? m->rows : m->cols;

    for (int i = 0; i < limit; i++) {
        int64_t val = m->matrix[i][i];
        if (val == 0)
            continue;

        rank++;

        int64_t abs_val = i64_abs(val);
        if (abs_val > 1) {
            if (torsions != NULL && nt < max_tors) {
                /* Cap at INT_MAX for the int out-buffer used by callers. */
                torsions[nt] = (abs_val > (int64_t)2147483647)
                                   ? 2147483647
                                   : (int)abs_val;
            }
            nt++;
        }
    }

    if (free_rank) *free_rank = rank;
    if (n_tors) *n_tors = nt;
    return 0;
}

int toSmithForm_extract(Mat* m,
                        int* free_rank,
                        int* torsions,
                        int max_tors,
                        int* n_tors)
{
    if (!m) return -1;
    toSmithForm(m);
    return smith_extract(m, free_rank, torsions, max_tors, n_tors);
}
