#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
typedef struct IntMatrix Mat;
struct IntMatrix {
  int rows;
  int cols;
  int** matrix; //2d pointer
  Mat* prev; //pointer to the previous matrix in the chain complex
  Mat* next; //next
  int* source;
  int* target;
};
Mat* createMat(int rows, int cols) {
  Mat* m = (Mat*)malloc(sizeof(Mat)); //storage management
  if (!m) return NULL;
  m->rows = rows;
  m->cols = cols;
  m->matrix = (int**)malloc(rows*sizeof(int*));
  //initialize matrix
  for (int i=0;i<rows;i++){
    m->matrix[i] = (int*)calloc(cols, sizeof(int));
  }
  m->prev = NULL;
  m->next = NULL;
  m->source = NULL;
  m->target = NULL;
  return m;
}
bool isDiag(Mat*m){
  //checks for diagonal matrix
  int rows = m->rows, cols = m->cols;
  for (int i=0;i<rows;i++){
    for (int j=0;j<cols;j++){
      //if any off-diagonal element is not a 0, its not a diagonal matrix
      if (i!=j && m->matrix[i][j]!=0){
        return false;
      }
    }
  }
  return true;
}
bool isnull(Mat* m){
  //all entries must be 0, so we stop if that isn't the case
  int rows = m->rows, cols = m->cols;
  for (int i=0;i<rows;i++){
    for (int j=0;j<cols;j++){
      if (m->matrix[i][j]!=0){
        return false;
      }
    }
  }
  return true;
}
//basic swapping
static void swaprows(Mat*m, int row1, int row2){
  int* temp = m->matrix[row1];
  m->matrix[row1] = m->matrix[row2];
  m->matrix[row2] = temp;
}
static void swapcols(Mat* m, int col1, int col2){
  for (int i=0;i<m->rows;i++){
    int temp =m->matrix[i][col1];
    m->matrix[i][col1] = m->matrix[i][col2];
    m->matrix[i][col2] = temp;
  }
}
//Consider our chain of complexes: When we swap the rows of one differential matrix d_n,
//since d_n is premultiplied onto the column of formally graded smoothings,
//we end up swapping the rows of the resultant output vector C_{n+1} in the same way we did for d_n
//Consider now the product d_{n+1}C_{n+1}, which will be different since every single dot product between the rows of
//d_{n+1} and the columns of C_{n+1} will change (that is, the one single column of C_{n+1} will change). 
//Thus, we need to permute the columns of d_{n+1} as well
//and this is indeed what takes place in the form of our second SwapRows.
//C_n --(d_n)--> C_{n+1} --(d_{n+1})--> C_{n+2}
void SwapRows(Mat* m, int row1, int row2){
  swaprows(m, row1, row2);
  if (m->next != NULL){
    swapcols(m->next, row1, row2);
  }
  if (m->target != NULL){
    int temp = m->target[row1];
    m->target[row1] = m->target[row2];
    m->target[row2] = temp;
  }
}
void SwapCols(Mat* m, int col1, int col2){
  swapcols(m, col1, col2);
  if (m->prev != NULL){
    swaprows(m->prev, col1, col2);
  }
  if (m->source != NULL){
    int temp = m->source[col1];
    m->source[col1] = m->source[col2];
    m->source[col2] = temp;
  }
}
static void addRow2(Mat* m, int a, int b, int n){
  for (int i=0;i<m->cols;i++){
    m->matrix[a][i] += (m->matrix[b][i]*n);
  }
}
static void addColumn2(Mat* m, int a, int b, int n) { 
  for (int i=0;i<m->rows;i++) {
    m->matrix[i][a] += (m->matrix[i][b]*n);
  }
}
static void multRow2(Mat* m, int a, int n) { 
  for (int i=0;i<m->cols;i++) {
      m->matrix[a][i]*=n;
  }
}
static void multColumn2(Mat* m, int a, int n) { 
  for (int i=0;i<m->rows;i++) {
      m->matrix[i][a]*=n;
  }
}
void addRow(Mat* m, int a, int b, int n){
  addRow2(m, a, b, n);
}
void addColumn(Mat* m, int a, int b, int n){
  addColumn2(m, a, b, n);
}
void multRow(Mat* m, int a, int n){
  multRow2(m, a, n);
}
void multColumn(Mat* m, int a, int n){
  multColumn2(m, a, n);
}
//query functions
int rowNonZeroes(Mat* m, int row){
  int count=0;
  for (int j=0;j<m->cols;j++){
    if (m->matrix[row][j]!=0){
      count++;
    }
  }
  return count;
}
int columnNonZeroes(Mat* m, int col){
  int count=0;
  for (int i=0;i<m->rows;i++){
    if (m->matrix[i][col]!=0){
      count++;
    }
  }
  return count;
}
int zeroRowsToEnd(Mat* m){
  int nzrows = m->rows;
  for (int i=0;i<nzrows;i++){
    while (rowNonZeroes(m, i) == 0 && i<nzrows){
      swaprows(m, i, nzrows-1);
      nzrows--;
    }
  }return nzrows;
}
int zeroColumnsToEnd(Mat* m){
  int nzcols = m->cols;
  for (int i=0;i<nzcols;i++){
    while (columnNonZeroes(m, i) == 0 && i<nzcols){
      swapcols(m, i, nzcols-1);
      nzcols--;
    }
  }return nzcols;
}
void toSF(Mat* m){
  for (int row=0,col=0; row<m->rows && col<m->cols; row++, col++) {
        while (row<m->rows && rowNonZeroes(m, row) == 0) row++;
        while (col<m->cols && columnNonZeroes(m, col) == 0) col++;
        if (row >= m->rows || col >= m->cols) break;
        if (row>col) {
            SwapRows(m, row, col);
            row = col;
        }else if (col>row) {
            SwapCols(m, row, col);
            col = row;
        }
        while (rowNonZeroes(m, row) != 1 || columnNonZeroes(m, col)!= 1 || m->matrix[row][col] <= 0) {
            for (int j=row; j<m->rows; j++) {
                if (m->matrix[j][col] < 0) {
                    multRow(m, j, -1);
                }
            }
            //euclidean algorithm on the column
            while (columnNonZeroes(m,col) != 1 || m->matrix[row][col]==0) {
                int min = -1;
                int idxmin = -1;
                for (int j=row; j<m->rows; j++) {
                    if ((m->matrix[j][col] < min || min == -1) && m->matrix[j][col] > 0) {
                        min = m->matrix[j][col];
                        idxmin = j;
                    }
                }
                if (idxmin != row) {
                    SwapRows(m,row, idxmin);
                }
                for (int j=row+1; j<m->rows; j++) {
                    if (m->matrix[j][col] != 0) {
                        addRow(m, j, row, -(m->matrix[j][col] / min));
                    }
                }
            }
            for (int j=col+1; j<m->cols; j++) {
                if (m->matrix[row][j] < 0) {
                    multColumn(m, j, -1);
                }
            }
            while (rowNonZeroes(m, row) != 1 || m->matrix[row][col] == 0) {
                int min = -1;
                int idxmin = -1;
                for (int j=col; j<m->cols; j++) {
                    if ((m->matrix[row][j]<min || min==-1) && m->matrix[row][j] > 0) {
                        min = m->matrix[row][j];
                        idxmin = j;
                    }
                }
                if (idxmin != col) {
                    SwapCols(m, col, idxmin);
                }
                for (int j=col+1; j<m->cols; j++) {
                    if (m->matrix[row][j] != 0) {
                        addColumn(m, j, col, -(m->matrix[row][j]/ min));
                    }
                }
            }
        }
    }
    zeroRowsToEnd(m);
    zeroColumnsToEnd(m);
}