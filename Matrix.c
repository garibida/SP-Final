#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define True 1
#define False 0

#define MatrixIterRows(A, i) for ((i) = 0; (i) < ((A) -> rows); (i)++)
#define MatrixIterCols(A, j) for ((j) = 0; (j) < ((A) -> cols); (j)++)
#define MatrixIterColsSym(A, i, j) for ((j) = 0; (j) <= (i); (j)++)

typedef double** Matrix_data;
typedef int bool;

typedef struct 
{
    Matrix_data data;
    int rows;
    int cols;
    bool isSymmetric;
} Matrix;

void freeMatrix(Matrix* A);
Matrix* createMatrix(int rows, int cols, bool isSymmetric);
Matrix_data createMatrixData(int rows, int cols);
Matrix_data createSymmetricMatrixData(int rows);
void updateMatrixSymmertircStatus(Matrix* A);
double getMatrixValue(Matrix* A, int row, int col);
void setMatrixValue(Matrix* A, int row, int col, double value);
void multiply_scalar(Matrix *A, double scalar);
Matrix* add(Matrix *A, Matrix *B);
Matrix* sub(Matrix *A, Matrix *B);
Matrix* multiply(Matrix* A, Matrix* B);
void printMatrix(Matrix* A);

void freeMatrix(Matrix* A) {
    int i;
    for (i = 0; i < A -> rows; i++) {
        free((A -> data)[i]);
    }
    free(A -> data);
    free(A);
}

Matrix* createMatrix(int rows, int cols, bool isSymmetric) {
    Matrix* A = (Matrix*) malloc(sizeof(Matrix));
    assert(A != NULL);
    A -> rows = rows;
    A -> cols = cols;
    A -> isSymmetric = isSymmetric;
    A -> data = isSymmetric ? createSymmetricMatrixData(rows) : createMatrixData(rows, cols);
    return A;
}

Matrix_data createMatrixData(int rows, int cols) {
    Matrix_data data;
    int i;

    data = (Matrix_data) calloc(rows, sizeof(double*));
    assert(data != NULL);
    for (i = 0; i < rows; i++) {
        data[i] = (double*) calloc(cols, sizeof(double));
        assert(data[i] != NULL);
    }

    return data;
}

Matrix_data createSymmetricMatrixData(int rows) {
    Matrix_data data;
    int i;

    data = (Matrix_data) calloc(rows, sizeof(double*));
    assert(data != NULL);
    for (i = 0; i < rows; i++) {
        data[i] = (double*) calloc(rows, sizeof(double));
        assert(data[i] != NULL);
    }

    return data;
}

void updateMatrixSymmertircStatus(Matrix* A) {
    Matrix_data data = A -> data;
    int i, j;
    bool isSym = True;
    MatrixIterRows(A, i) {
        MatrixIterColsSym(A, i, j) {
            if(data[i][j] != data[j][i]) {
                isSym = False;
            }
        }
    }
    A -> isSymmetric = isSym;
}

double getMatrixValue(Matrix* A, int row, int col) {
    assert(row < A -> rows && col < A -> cols);
    return (A -> isSymmetric && col > row) ?
                (A -> data)[col][row]:
                (A -> data)[row][col];
}

void setMatrixValue(Matrix* A, int row, int col, double value) {
    assert(row < A -> rows && col < A -> cols);
    if (A -> isSymmetric && col > row) {
        (A -> data)[col][row] = value;
    } else {
        (A -> data)[row][col] = value;
    }
}

void multiply_scalar(Matrix *A, double scalar) {
    Matrix_data data = A -> data;
    int i, j;
    MatrixIterRows(A, i) {
        if (A -> isSymmetric){
            MatrixIterColsSym(A, i, j) {
                data[i][j] = data[i][j] * scalar;
            }
        }
        else {
            MatrixIterCols(A, j) {
                data[i][j] = data[i][j] * scalar;
            }
        }
    }
}

Matrix* add(Matrix *A, Matrix *B) {
    Matrix* C;
    int i, j;
    bool isSymmetric;
    assert(A -> rows == B -> rows);
    assert(A -> cols == B -> cols);

    isSymmetric = A -> isSymmetric && B -> isSymmetric;

    C = createMatrix(A -> rows, B -> cols, isSymmetric);

    MatrixIterRows(A, i) {
        if (C -> isSymmetric){
            MatrixIterColsSym(A, i, j) {
                (C -> data)[i][j] = getMatrixValue(A,i,j) + getMatrixValue(B,i,j);
            }
        }
        else {
            MatrixIterCols(A, j) {
                (C -> data)[i][j] = getMatrixValue(A,i,j) + getMatrixValue(B,i,j);
            }
        }
    }

    return C;
}

Matrix* sub(Matrix *A, Matrix *B) {
    Matrix* C;
    int i, j;
    bool isSymmetric;
    assert(A -> rows == B -> rows);
    assert(A -> cols == B -> cols);

    isSymmetric = A -> isSymmetric && B -> isSymmetric;

    C = createMatrix(A -> rows, B -> cols, isSymmetric);

    MatrixIterRows(A, i) {
        if (C -> isSymmetric){
            MatrixIterColsSym(A, i, j) {
                (C -> data)[i][j] = getMatrixValue(A,i,j) - getMatrixValue(B,i,j);
            }
        }
        else {
            MatrixIterCols(A, j) {
                (C -> data)[i][j] = getMatrixValue(A,i,j) - getMatrixValue(B,i,j);
            }
        }
    }

    return C;
}

Matrix* multiply(Matrix* A, Matrix* B) {
    Matrix* C;
    int i, j, k;
    double value;
    assert(A -> cols == B -> rows);

    C = createMatrix(A -> rows, B -> cols, False);
    MatrixIterRows(C, i) {
        MatrixIterCols(C, j) {
            value = 0;
            MatrixIterCols(A, k) {
                value += getMatrixValue(A,i,k) * getMatrixValue(B,k,j);
            }
            (C -> data)[i][j] = value;
        }
    }

    return C;
}

void printMatrix(Matrix* A) {
    int i, j;
    printf("===================\n");
    for (i = 0; i < A -> rows; i++)
    {
        for(j = 0; j< A -> cols; j++)
        {
            printf("%f     ", getMatrixValue(A,i,j));
        }
        printf("\n");
    }
}

int main() {
    Matrix* B, *C, *A;
    A = createMatrix(3, 3, True);
    setMatrixValue(A, 0 ,0, 1.0);
    setMatrixValue(A, 1 ,0, 2.0);
    setMatrixValue(A, 1 ,1, 3.0);
    setMatrixValue(A, 2 ,0, 4.0);
    setMatrixValue(A, 2 , 1, 5.0);
    setMatrixValue(A, 2 , 2, 6.0);
    
    B = createMatrix(3, 2, False);
    (B -> data) [0][0] = 1.0;
    (B -> data) [1][0] = 8.0;
    (B -> data) [2][0] = 7.0;
    (B -> data) [0][1] = 11.0;
    (B -> data) [1][1] = 6.0;
    (B -> data) [2][1] = 1.0;

    printMatrix(A);
    printMatrix(B);

    C = multiply(A, B);

    printMatrix(C);

    freeMatrix(A);
    freeMatrix(B);
    freeMatrix(C);

    return 0;
}