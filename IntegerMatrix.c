#include "IntegerMatrix.h"
#include <stdio.h>
/* * Allocates memory for the matrix structure and its internal 2D grid.
 * Initializes tracking pointers and links to NULL.
 */
Mat* createMat(int rows, int cols) {
    Mat* m = (Mat*)malloc(sizeof(Mat));
    if (!m) return NULL;
    
    m->rows = rows;
    m->cols = cols;
    
    // Allocate array of pointers for rows using 64-bit type
    m->matrix = (int64_t**)malloc(rows * sizeof(int64_t*));
    for (int i = 0; i < rows; i++) {
        // Allocate and zero-initialize columns for each row
        m->matrix[i] = (int64_t*)calloc(cols, sizeof(int64_t));
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
    
    // Free each individual row row-by-row
    for (int i = 0; i < m->rows; i++) {
        free(m->matrix[i]);
    }
    free(m->matrix);
    
    // Free generator tracking arrays if they were allocated
    if (m->source != NULL) free(m->source);
    if (m->target != NULL) free(m->target);
    
    free(m);
}

/* * Isolates the matrix by breaking its links within the chain complex.
 */
void isolateMat(Mat* m) {
    if (!m) return;
    
    // Safely disconnect from the neighbors in the chain complex
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
        m->matrix[a][i] += (m->matrix[b][i] * n);
    }
}

/* Low-level operation: Column Addition -> Col(a) = Col(a) + n * Col(b) */
void addColumn2(Mat* m, int a, int b, int64_t n) {
    for (int i = 0; i < m->rows; i++) {
        m->matrix[i][a] += (m->matrix[i][b] * n);
    }
}

/* Low-level operation: Multiplies a row by a scalar factor */
void multRow2(Mat* m, int a, int64_t n) {
    for (int i = 0; i < m->cols; i++) {
        m->matrix[a][i] *= n;
    }
}

/* Low-level operation: Multiplies a column by a scalar factor */
void multColumn2(Mat* m, int a, int64_t n) {
    for (int i = 0; i < m->rows; i++) {
        m->matrix[i][a] *= n;
    }
}

/* * High-level operation: Cascading Row Swap.
 * Changing the basis of the target space of d_n (swapping rows) requires 
 * swapping the columns of d_{n+1} to maintain the chain complex property (d ◦ d = 0).
 */
void SwapRows(Mat* m, int row1, int row2) {
    if (row1 == row2) return;
    swaprows(m, row1, row2);
    
    // Cascade to the next operator's columns to maintain algebraic consistency
    if (m->next != NULL) {
        swapcols(m->next, row1, row2);
    }
    
    // Update target generator indices mappings
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
    
    // Cascade to the previous operator's rows
    if (m->prev != NULL) {
        swaprows(m->prev, col1, col2);
    }
    
    // Update source generator indices mappings
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
    // Isolate matrix to prevent cascading actions from breaking the rest of the chain complex
    isolateMat(m);

    int row = 0, col = 0;
    while (row < m->rows && col < m->cols) {
        // Find next active rows and columns that are non-zero
        while (row < m->rows && rowNonZeroes(m, row) == 0) row++;
        while (col < m->cols && columnNonZeroes(m, col) == 0) col++;
        
        if (row >= m->rows || col >= m->cols) break;
        
        // Bring pivot element to matching diagonal position
        if (row > col) { SwapRows(m, row, col); row = col; }
        else if (col > row) { SwapCols(m, row, col); col = row; }
        
        // Euclidean elimination loop for column clearing
        while (rowNonZeroes(m, row) != 1 || columnNonZeroes(m, col) != 1 || m->matrix[row][col] <= 0) {
            // Force values in current column to be positive
            for (int j = row; j < m->rows; j++) {
                if (m->matrix[j][col] < 0) {
                    multRow(m, j, -1);
                }
            }
            
            // Eliminate all entries under the pivot inside the column using the Euclidean algorithm
            while (columnNonZeroes(m, col) != 1 || m->matrix[row][col] == 0) {
                int64_t min_val = -1;
                int idxmin = -1;
                
                // Locate the absolute minimum positive element to serve as the new temporary pivot
                for (int j = row; j < m->rows; j++) {
                    int64_t v = m->matrix[j][col];
                    if (v > 0 && (min_val == -1 || v < min_val)) {
                        min_val = v; idxmin = j;
                    }
                }
                
                if (idxmin != row) SwapRows(m, row, idxmin);
                
                // Reduce other rows using integer division (Euclidean Step)
                for (int j = row + 1; j < m->rows; j++) {
                    if (m->matrix[j][col] != 0) {
                        addRow(m, j, row, -(m->matrix[j][col] / min_val));
                    }
                }
            }
        }
        
        // Force values in the row to be positive
        for (int j = col + 1; j < m->cols; j++) {
            if (m->matrix[row][j] < 0) {
                multColumn(m, j, -1);
            }
        }
        
        // Eliminate entries to the right of the pivot inside the row using the Euclidean algorithm
        while (rowNonZeroes(m, row) != 1 || m->matrix[row][col] == 0) {
            int64_t min_val = -1;
            int idxmin = -1;
            
            for (int j = col; j < m->cols; j++) {
                int64_t v = m->matrix[row][j];
                if (v > 0 && (min_val == -1 || v < min_val)) {
                    min_val = v; idxmin = j;
                }
            }
            
            if (idxmin != col) SwapCols(m, col, idxmin);
            
            // Reduce remaining columns
            for (int j = col + 1; j < m->cols; j++) {
                if (m->matrix[row][j] != 0) {
                    addColumn(m, j, col, -(m->matrix[row][j] / min_val));
                }
            }
        }
    }
    
    // Clean up structure by packing all zero dimensions together
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