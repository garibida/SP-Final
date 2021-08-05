#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define True 1
#define False 0

#define MatrixIterRows(A, i) for ((i) = 0; (i) < ((A) -> rows); (i)++)
#define MatrixIterCols(A, j) for ((j) = 0; (j) < ((A) -> cols); (j)++)


typedef double** Matrix_data;
typedef int bool;

typedef struct 
{
    Matrix_data data;
    int rows;
    int cols;
} Matrix;

void freeMatrix(Matrix* A) {
    int i;
    for (i = 0; i < A -> rows; i++) {
        free((A -> data)[i]);
    }
    free(A -> data);
    free(A);
}

Matrix* createMatrix(int rows, int cols) {
    Matrix* A;
    Matrix_data data;
    int i;

    A = (Matrix*) malloc(sizeof(Matrix));
    assert(A != NULL);
    A -> rows = rows;
    A -> cols = cols;
    data = (Matrix_data) calloc(rows, sizeof(double*));
    assert(data != NULL);
    for (i = 0; i < rows; i++) {
        data[i] = (double*) calloc(cols, sizeof(double));
        assert(data[i] != NULL);
    }
    A -> data = data;

    return A;
}

void multiply_scalar(Matrix *A, double scalar) {
    Matrix_data data = A -> data;
    int i, j;
    MatrixIterRows(A, i) {
        MatrixIterCols(A, j) {
            data[i][j] = data[i][j] * scalar;
        }
    }
}

Matrix* add(Matrix *A, Matrix *B) {
    Matrix* C;
    int i, j;
    assert(A -> rows == B -> rows);
    assert(A -> cols == B -> cols);

    C = createMatrix(A -> rows, B -> cols);

    MatrixIterRows(A, i) {
        MatrixIterCols(A, j) {
            (C -> data)[i][j] = (A -> data)[i][j] + (B -> data)[i][j];
        }
    }

    return C;
}

Matrix* sub(Matrix *A, Matrix *B) {
    Matrix* C;
    int i, j;
    assert(A -> rows == B -> rows);
    assert(A -> cols == B -> cols);

    C = createMatrix(A -> rows, B -> cols);

    MatrixIterRows(A, i) {
        MatrixIterCols(A, j) {
            (C -> data)[i][j] = (A -> data)[i][j] - (B -> data)[i][j];
        }
    }

    return C;
}

/*
Matrix* multiply(Matrix A, Matrix B) {
    
}
*/

void printMatrix(Matrix* A) {
    int i, j;
    for (i = 0; i < A -> rows; i++)
    {
        for(j = 0; j< A -> cols; j++)
        {
            printf("%f     ", (A -> data)[i][j]);
        }
        printf("\n");
    }
}

int main() {
    /*
    Matrix* B, *C, *A = createMatrix(2, 2);
    (A -> data) [0][0] = 1.0;
    (A -> data) [1][0] = 2.0;
    (A -> data) [0][1] = 3.0;
    (A -> data) [1][1] = 4.0;

    
    B = createMatrix(2, 2);
    (B -> data) [0][0] = 1.0;
    (B -> data) [1][0] = 0;
    (B -> data) [0][1] = 0;
    (B -> data) [1][1] = 1.0;

    printMatrix(A);
    printMatrix(B);

    C = add(A, B);

    printMatrix(C);

    freeMatrix(A);
    freeMatrix(B);
    freeMatrix(C);
    */

    return 0;
}