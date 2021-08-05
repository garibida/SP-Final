#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h> /* check */

#define true 1
#define false 0

#define MatrixIterRows(A, i) for ((i) = 0; (i) < ((A) -> rows); (i)++)
#define MatrixIterCols(A, j) for ((j) = 0; (j) < ((A) -> cols); (j)++)
#define MatrixIterColsSym(A, i, j) for ((j) = 0; (j) <= (i); (j)++)

typedef double** Matrix_data;
typedef double* Point_data; 
typedef int bool;

typedef struct {
    Matrix_data data;
    int rows;
    int cols;
    bool isSymmetric;
} Matrix;

typedef struct {
    Point_data data;
    int d;
} Point;


/* ! validated, ? coded */ 
int validateInput(); /* !? */ 

/* Algorith's operations section */
Matrix* computeMatrixW(Point** pointsArr, int n);
Matrix* computeMatrixD(Point** pointsArr); /* !? */ 
Matrix* computeMatrixL(Point** pointsArr); /* !? */ 
Matrix* computeMatrixLnorm(Point** pointsArr); /* !? */ 

/* Point's operations section */
Point* createPoint(int d);
void printPoint(Point* point);
void printPointsArr(Point **pointArr, int n);
int isPointsEquel(Point* point1, Point* point2);
double computeDist(Point* point1, Point* point2);
double computeDistW(Point* point1, Point* point2); 
void freeMemPoint(Point* point);
void freeMemPointsArr(Point **pointsArr, int n);
/* double* copy_point(double* point, int d); */ 

/* Matrix's operations section */
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
bool isMatrixEqual(Matrix *A, Matrix *B);
void printMatrix(Matrix* A);

/* test section */
void testMain(bool isdebug);
void testCalcMatrixW(bool isdebug);
void testMultiplyMatrixs();

/* ######################### */
/* Matrix operations section */
/* ######################### */

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
    MatrixIterRows(A, i) {
        MatrixIterColsSym(A, i, j) {
            if(data[i][j] != data[j][i]) {
                A -> isSymmetric = false;
                return;
            }
        }
    }
    A -> isSymmetric = true;
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

    C = createMatrix(A -> rows, B -> cols, false);
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

void freeMatrix(Matrix* A) {
    int i;
    for (i = 0; i < A -> rows; i++) {
        free((A -> data)[i]);
    }
    free(A -> data);
    free(A);
}

bool isMatrixEqual(Matrix *A, Matrix *B) {
    int i, j;
    double epsilon = 0.00000001;
    bool isSymmetric;
    assert(A -> rows == B -> rows);
    assert(A -> cols == B -> cols);

    isSymmetric = A -> isSymmetric && B -> isSymmetric;

    MatrixIterRows(A, i) {
        if (isSymmetric){
            MatrixIterColsSym(A, i, j) {
                if(fabs(getMatrixValue(A,i,j) - getMatrixValue(B,i,j)) > epsilon){
                    return false;
                }
            }
        }
        else {
            MatrixIterCols(A, j) {
                if(fabs(getMatrixValue(A,i,j) - getMatrixValue(B,i,j)) > epsilon){
                    return false;
                }
            }
        }
    }

    return true;
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

/* ######################## */
/* Points operation section */
/* ######################## */

Point* createPoint(int d) {
    Point *point;
    Point_data data;
    
    point = (Point*) malloc(sizeof(Point));
    assert(point != NULL);
    point -> d = d;

    data = (Point_data) calloc(d, sizeof(double));
    assert(data != NULL);
    point -> data = data;

    return point;
}

void printPoint(Point *point) {
    int i, dim;
    dim = point -> d;
    for (i = 0; i < dim; i++) {
        printf("%.4f", point -> data[i]);
        if (i != dim - 1) {
            printf(",");
        }
    }
    printf("\n");
}

void printPointsArr(Point **pointArr, int n) {
    int i;
    for (i = 0; i < n; i++) {
        printPoint(pointArr[i]);
    }
}

int isPointsEquel(Point *point1, Point *point2){
    int i;
    for (i = 0; i < point1 -> d; i++) {
        if (point1 -> data[i] != point2 -> data[i]) {
            return false;
        }
    }
    return true;
} 

void freeMemPoint(Point *point) {
    if (point != NULL) {
        free(point -> data);
        free(point); 
    }
}

void freeMemPointsArr(Point **pointsArr, int n) {
    int i; 
    for (i = 0; i < n; i++) {
        freeMemPoint(pointsArr[i]);
    }
    free(pointsArr);
}

/* ######### */
/* Algorithm */
/* ######### */

double computeDist(Point *point1, Point *point2) {
    double dist = 0, tmp = 0; 
    int i;
    for (i = 0; i < point1 -> d; i++) {
        tmp = (point1 -> data[i] - point2 -> data[i]);
        dist += tmp * tmp;
    }
    return sqrt(dist);
}

double computeDistW(Point *point1, Point *point2) {
    double dist = computeDist(point1, point2);
    return exp(-0.5 * dist * dist);
}

Matrix* computeMatrixW(Point **pointsArr, int n) {
    Matrix *W; 
    Point *pointI, *pointJ;
    int i, j;
    double wTmp;
    W = createMatrix(n,n, true); 
    
    MatrixIterRows(W, i) {
        pointI = pointsArr[i]; 
        MatrixIterColsSym(W, i, j) {
            pointJ = pointsArr[j];
            if (i == j) {
                W -> data[i][j] = 0;
            } else {
                wTmp = computeDistW(pointI, pointJ);
                W -> data[i][j] = wTmp;
                W -> data[j][i] = wTmp;
            }
        }
    }

    return W;
}

/* ############# */
/* Tests section */
/* ############# */

void testMain(bool isdebug) {
    testCalcMatrixW(isdebug);
    testMultiplyMatrixs();
}

void testCalcMatrixW(bool isdebug) {
    Matrix *W, *A;
    Point **pointsArr;
    int i, b1;
    pointsArr = calloc(3, sizeof(Point));
    assert(pointsArr != NULL); 
    pointsArr[0] = createPoint(3); 
    pointsArr[1] = createPoint(3); 
    for (i = 0; i < 3; i++) {
        pointsArr[1] -> data[i] = 1.0;
    }
    
    pointsArr[2] = createPoint(3); 
    for (i = 0; i < 3; i++) {
        pointsArr[2] -> data[i] = 2.0;
    }
    
    A = createMatrix(3, 3, true);
    setMatrixValue(A, 0 ,0, 0.0);
    setMatrixValue(A, 0 ,1, exp(-1.5));
    setMatrixValue(A, 0 ,2, exp(-6));
    setMatrixValue(A, 1 ,0, exp(-1.5));
    setMatrixValue(A, 1 ,1, 0.0);
    setMatrixValue(A, 1 ,2, exp(-1.5));
    setMatrixValue(A, 2 ,0, exp(-6));
    setMatrixValue(A, 2 ,1, exp(-1.5));
    setMatrixValue(A, 2 ,2, 0.0);

    W = computeMatrixW(pointsArr, 3);
    
    if (isdebug == 1) {
        printf("points array: \n");
        printPointsArr(pointsArr, 3);
        printf("Matrix W calculated: \n");
        printMatrix(W);
        printf("Matrix A correct Matrix\n");
        printMatrix(A);
    }
    
    b1 = isMatrixEqual(W,A);
    printf("test name: 'testCalcMatrixW'\tresult: %d\n", b1);

    freeMatrix(W);
    freeMatrix(A);
    freeMemPointsArr(pointsArr, 3);
}

void testMultiplyMatrixs() {
    Matrix* B, *C, *A;
    A = createMatrix(3, 3, true);
    setMatrixValue(A, 0 ,0, 1.0);
    setMatrixValue(A, 1 ,0, 2.0);
    setMatrixValue(A, 1 ,1, 3.0);
    setMatrixValue(A, 2 ,0, 4.0);
    setMatrixValue(A, 2 , 1, 5.0);
    setMatrixValue(A, 2 , 2, 6.0);
    
    B = createMatrix(3, 2, false);
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
}

int main() {
    testMain(true);
    return 0;
}