#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
typedef struct IntMatrix Mat;
struct IntMatrix {
    int rows;
    int cols;
    int** matrix;
    Mat* prev;
    Mat* next;
    int* source;
    int* target;
};
Mat* createMat(int rows, int cols) {
    Mat* m = (Mat*)malloc(sizeof(Mat));
    if (!m) return NULL;
    m->rows = rows;
    m->cols = cols;
    m->matrix = (int**)malloc(rows * sizeof(int*));
    for (int i = 0; i < rows; i++) {
        m->matrix[i] = (int*)calloc(cols, sizeof(int));
    }
    m->prev = NULL;
    m->next = NULL;
    m->source = NULL;
    m->target = NULL;
    return m;
}
void freeMat(Mat* m) {
    if (!m) return;
    for (int i = 0; i < m->rows; i++) free(m->matrix[i]);
    free(m->matrix);
    free(m);
}
bool isDiag(Mat* m) {
    int rows = m->rows, cols = m->cols;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (i != j && m->matrix[i][j] != 0)
                return false;
    return true;
}
bool isnull(Mat* m) {
    int rows = m->rows, cols = m->cols;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (m->matrix[i][j] != 0)
                return false;
    return true;
}
int rowNonZeroes(Mat* m, int row) {
    int count = 0;
    for (int j = 0; j < m->cols; j++)
        if (m->matrix[row][j] != 0) count++;
    return count;
}
int columnNonZeroes(Mat* m, int col) {
    int count = 0;
    for (int i = 0; i < m->rows; i++)
        if (m->matrix[i][col] != 0) count++;
    return count;
}
void swaprows(Mat* m, int row1, int row2) {
    int* temp = m->matrix[row1];
    m->matrix[row1] = m->matrix[row2];
    m->matrix[row2] = temp;
}
void swapcols(Mat* m, int col1, int col2) {
    for (int i = 0; i < m->rows; i++) {
        int temp = m->matrix[i][col1];
        m->matrix[i][col1] = m->matrix[i][col2];
        m->matrix[i][col2] = temp;
    }
}
void addRow2(Mat* m, int a, int b, int n) {
    for (int i = 0; i < m->cols; i++)
        m->matrix[a][i] += (m->matrix[b][i] * n);
}
void addColumn2(Mat* m, int a, int b, int n) {
    for (int i = 0; i < m->rows; i++)
        m->matrix[i][a] += (m->matrix[i][b] * n);
}
void multRow2(Mat* m, int a, int n) {
    for (int i = 0; i < m->cols; i++)
        m->matrix[a][i] *= n;
}
void multColumn2(Mat* m, int a, int n) {
    for (int i = 0; i < m->rows; i++)
        m->matrix[i][a] *= n;
}
void SwapRows(Mat* m, int row1, int row2) {
    if (row1 == row2) return;
    swaprows(m, row1, row2);
    if (m->next != NULL)
        swapcols(m->next, row1, row2);
    if (m->target != NULL) {
        int temp = m->target[row1];
        m->target[row1] = m->target[row2];
        m->target[row2] = temp;
    }
}
void SwapCols(Mat* m, int col1, int col2) {
    if (col1 == col2) return;
    swapcols(m, col1, col2);
    if (m->prev != NULL)
        swaprows(m->prev, col1, col2);
    if (m->source != NULL) {
        int temp = m->source[col1];
        m->source[col1] = m->source[col2];
        m->source[col2] = temp;
    }
}
void addRow(Mat* m, int a, int b, int n) { addRow2(m, a, b, n); }
void addColumn(Mat* m, int a, int b, int n) { addColumn2(m, a, b, n); }
void multRow(Mat* m, int a, int n) { multRow2(m, a, n); }
void multColumn(Mat* m, int a, int n) { multColumn2(m, a, n); }
int zeroRowsToEnd(Mat* m) {
    int nzrows = m->rows;
    for (int i = 0; i < nzrows; i++)
        while (rowNonZeroes(m, i) == 0 && i < nzrows)
            SwapRows(m, i, --nzrows);
    return nzrows;
}
int zeroColumnsToEnd(Mat* m) {
    int nzcols = m->cols;
    for (int i = 0; i < nzcols; i++)
        while (columnNonZeroes(m, i) == 0 && i < nzcols)
            SwapCols(m, i, --nzcols);
    return nzcols;
}
void toSmithForm(Mat* m) {
    int row = 0, col = 0;
    while (row < m->rows && col < m->cols) {
        while (row < m->rows && rowNonZeroes(m, row) == 0) row++;
        while (col < m->cols && columnNonZeroes(m, col) == 0) col++;
        if (row >= m->rows || col >= m->cols) break;
        if (row > col) { SwapRows(m, row, col); row = col; }
        else if (col > row) { SwapCols(m, row, col); col = row; }
        while (rowNonZeroes(m, row) != 1 || columnNonZeroes(m, col) != 1 || m->matrix[row][col] <= 0) {
            for (int j = row; j < m->rows; j++)
                if (m->matrix[j][col] < 0)
                    multRow(m, j, -1);
            while (columnNonZeroes(m, col) != 1 || m->matrix[row][col] == 0) {
                int min_val = -1, idxmin = -1;
                for (int j = row; j < m->rows; j++) {
                    int v = m->matrix[j][col];
                    if (v > 0 && (min_val == -1 || v < min_val)) {
                        min_val = v; idxmin = j;
                    }
                }
                if (idxmin != row) SwapRows(m, row, idxmin);
                for (int j = row + 1; j < m->rows; j++)
                    if (m->matrix[j][col] != 0)
                        addRow(m, j, row, -(m->matrix[j][col] / min_val));
            }
        }
        for (int j = col + 1; j < m->cols; j++)
            if (m->matrix[row][j] < 0)
                multColumn(m, j, -1);
        while (rowNonZeroes(m, row) != 1 || m->matrix[row][col] == 0) {
            int min_val = -1, idxmin = -1;
            for (int j = col; j < m->cols; j++) {
                int v = m->matrix[row][j];
                if (v > 0 && (min_val == -1 || v < min_val)) {
                    min_val = v; idxmin = j;
                }
            }
            if (idxmin != col) SwapCols(m, col, idxmin);
            for (int j = col + 1; j < m->cols; j++)
                if (m->matrix[row][j] != 0)
                    addColumn(m, j, col, -(m->matrix[row][j] / min_val));
        }
    }
    zeroRowsToEnd(m);
    zeroColumnsToEnd(m);
}
void printMat(Mat* m) {
    for (int i=0;i<m->rows;i++) {
        for (int j=0;j<m->cols;j++)
            printf("%4d",m->matrix[i][j]);
        printf("\n");
    }
}