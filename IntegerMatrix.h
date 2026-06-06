#ifndef INT_MATRIX_H
#define INT_MATRIX_H

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct Mat Mat;

struct Mat {
    int rows;
    int cols;
    int64_t** matrix; //replacement of Java's BigInteger[][]
    
    //chaining complexes
    Mat* prev;
    Mat* next;
    
    //tracking arrays
    //replacement for Java's List<Integer> source, target;
    int* source;
    int* target; 
};
//memory management
Mat* createMat(int rows, int cols);
void freeMat(Mat* m);
//gaussian elimination and smith
//column
void multColumn(Mat* m, int col, int64_t scalar);
void SwapCols(Mat* m, int col1, int col2);
void addColumn(Mat* m, int dest_col, int src_col, int64_t scalar);

// Row Operations (Derived from the Java IntMatrix behavior)
void multRow(Mat* m, int row, int64_t scalar);
void SwapRows(Mat* m, int row1, int row2);
void addRow(Mat* m, int dest_row, int src_row, int64_t scalar);

/* ========================================================================= */
/* Utility & Compaction Functions                                            */
/* ========================================================================= */

int rowNonZeroes(Mat* m, int row);
void zeroRowsToEnd(Mat* m);
void zeroColumnsToEnd(Mat* m);

#endif