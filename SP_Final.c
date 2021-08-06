#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h> /* check */

#define true 1
#define false 0
#define TestMode 1

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
Matrix* computeMatrixD(Matrix *W);
Matrix* computeMatrixDMinusHalf(Matrix *D);
Matrix* computeMatrixL(Matrix *W, Matrix *D); 
Matrix* computeMatrixLnorm(Matrix *L, Matrix *D); /* write tester */ 

/* Point's operations section */
Point* createPoint(int d);
Point* createPointWithVals(int d, double *values); /* !? */ 
Point* setDataPointVal(Point *point, int index, double value); /* !? */ 
double getDataPointVal(Point *point, int index); /* !? */ 
void printPoint(Point* point);
void printPointsArr(Point **pointArr, int n);
int isPointsEquel(Point *point1, Point* point2);
double computeDist(Point *point1, Point* point2);
double computeDistW(Point *point1, Point* point2); 
void freeMemPoint(Point *point);
void freeMemPointsArr(Point **pointsArr, int n);
Point* copy_point(Point *point);

/* Matrix's operations section */
void freeMatrix(Matrix* A);
Matrix* createMatrix(int rows, int cols, bool isSymmetric);
Matrix* createUnitMatrix(int n);
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
void testMain(bool isDebug);
Point** createPointsArrTest();
void testMultiplyMatrixs(bool isDebug);
void Test1(bool isDebug);
/*testMatrixLnorm(isDebug);*/

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

Matrix* createUnitMatrix(int dim) {
    Matrix* I; 
    int i;
    I = createMatrix(dim, dim, true);
    MatrixIterRows(I, i) {
        setMatrixValue(I, i, i, 1.0);
    }
    return I;
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
    int i, j;
    MatrixIterRows(A, i) {
        if (A -> isSymmetric){
            MatrixIterColsSym(A, i, j) {
                setMatrixValue(A, i, j, getMatrixValue(A,i,j) * scalar);
            }
        }
        else {
            MatrixIterCols(A, j) {
                setMatrixValue(A, i, j, getMatrixValue(A,i,j) * scalar);
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
            setMatrixValue(C, i, j, value);
        }
    }

    return C;
}

void freeMatrix(Matrix *A) {
    int i;
    MatrixIterRows(A, i) {
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
    MatrixIterRows(A,i) {
        MatrixIterCols(A,j) {
            printf("%.4f    ", getMatrixValue(A,i,j));
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

Point* createPointWithVals(int d, double *values) {
    /* add assertion for length comp to len(values) */
    Point *point;
    int i;
    
    point = createPoint(d);
    for (i = 0; i < d; i++) {
        setDataPointVal(point, i, values[i]);
    }
    return point;
}

Point* setDataPointVal(Point *point, int index, double value) {
    point->data[index] = value;
    return point;
}

double getDataPointVal(Point *point, int index) {
    return point->data[index];
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
        if (getDataPointVal(point1, i) != getDataPointVal(point2, i)) {
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

Point* copy_point(Point *point) {
    int i;
    double val;
    Point *newPoint = createPoint(point -> d); 
    for (i = 0; i < point -> d; i++) {
        val = getDataPointVal(point,i);
        setDataPointVal(newPoint, i, val);
    }
    return newPoint;
}

/* ######### */
/* Algorithm */
/* ######### */

double computeDist(Point *point1, Point *point2) {
    double dist = 0, tmp = 0; 
    int i;
    for (i = 0; i < point1 -> d; i++) {
        tmp = getDataPointVal(point1, i) - getDataPointVal(point2,i);
        dist += tmp * tmp;
    }
    return sqrt(dist);
}

double computeDistW(Point *point1, Point *point2) {
    double dist = computeDist(point1, point2);
    return exp(-0.5 * dist * dist);
}

Matrix* computeMatrixW(Point **pointsArr, int dim) {
    Matrix *W; 
    Point *pointI, *pointJ;
    int i, j;
    double wVal;
    W = createMatrix(dim, dim, true); 
    
    MatrixIterRows(W, i) {
        pointI = pointsArr[i]; 
        MatrixIterColsSym(W, i, j) {
            pointJ = pointsArr[j];
            if (i == j) {
                setMatrixValue(W, i, j, 0);
            } else {
                wVal = computeDistW(pointI, pointJ);
                setMatrixValue(W, i, j, wVal);
            }
        }
    }
    return W;
}

Matrix* computeMatrixD(Matrix *W) {
    Matrix *D; 
    int i, j, dim;
    double dVal = 0;
    dim = W -> cols;
    D = createMatrix(dim, dim, true); 
    
    MatrixIterRows(W, i) {
        MatrixIterCols(W, j) {
            dVal += getMatrixValue(W,i,j);
        }
        setMatrixValue(D,i,i,dVal);
        dVal = 0;
    }
    return D;
}

Matrix* computeMatrixDMinusHalf(Matrix *D) {
    Matrix *D2; 
    int i, dim;
    double dVal;
    dim = D -> cols;
    D2 = createMatrix(dim, dim, true); 
    
    MatrixIterRows(D, i) {
        dVal = 1 / (sqrt(getMatrixValue(D,i,i)));
        setMatrixValue(D2,i,i,dVal);
    }
    return D2;
}

Matrix* computeMatrixL(Matrix *W, Matrix *D) {
    return sub(D, W);
}

Matrix* computeMatrixLnorm(Matrix *W ,Matrix *D) {
    Matrix *D2, *I, *tmp;
    D2 = computeMatrixDMinusHalf(D);
    I = createUnitMatrix(W -> rows);
    tmp = multiply(multiply(D2, W), D2);
    return sub(I, tmp);
}

/* ############# */
/* Tests section */
/* ############# */

void testMain(bool isDebug) {
    testMultiplyMatrixs(isDebug);
    Test1(isDebug);
    /*testMatrixLnorm(isDebug);*/
}

void testMultiplyMatrixs(bool isDebug) {
    Matrix *A, *B, *C, *D;
    A = createMatrix(3, 3, true);
    setMatrixValue(A , 0 , 0 , 1.0);
    setMatrixValue(A , 1 , 0 , 2.0);
    setMatrixValue(A , 1 , 1 , 3.0);
    setMatrixValue(A , 2 , 0 , 4.0);
    setMatrixValue(A , 2 , 1 , 5.0);
    setMatrixValue(A , 2 , 2 , 6.0);
    
    B = createMatrix(3, 2, false);
    setMatrixValue(B , 0 , 0 , 1.0);
    setMatrixValue(B , 1 , 0 , 8.0);
    setMatrixValue(B , 2 , 0 , 7.0);
    setMatrixValue(B , 0 , 1 , 11.0);
    setMatrixValue(B , 1 , 1 , 6.0);
    setMatrixValue(B , 2 , 1 , 1.0);

    if (isDebug) {
        printf("Matrix A: \n");
        printMatrix(A);
        printf("Matrix B: \n");
        printMatrix(B);
    }
    
    C = multiply(A, B);
    if (isDebug) {
        printf("Matrix C calculted: \n");
        printMatrix(C);
    }

    D = createMatrix(3, 2, false);
    setMatrixValue(D , 0 , 0 , 45.0);
    setMatrixValue(D , 1 , 0 , 61.0);
    setMatrixValue(D , 2 , 0 , 86.0);
    setMatrixValue(D , 0 , 1 , 27.0);
    setMatrixValue(D , 1 , 1 , 45.0);
    setMatrixValue(D , 2 , 1 , 80.0);

    if (isDebug) {
        printf("Matrix D result: \n");
        printMatrix(D);
    }

    (isMatrixEqual(C,D)) ?
        printf("'test Multiply Matrixs'\t\tresult: Great!\n") : 
        printf("'test Multiply Matrixs'\t\tresult: Problem!\n");

    freeMatrix(A);
    freeMatrix(B);
    freeMatrix(C);
    freeMatrix(D);
}

Point** pointsForTest1() {
    Point **pointsArr;
    int numOfPoints = 3;
    double pointVal1[3] = {0,0,0};
    double pointVal2[3] = {1,1,1}; 
    double pointVal3[3] = {2,2,2};

    pointsArr = calloc(numOfPoints, sizeof(Point));
    assert(pointsArr != NULL); 
    
    pointsArr[0] = createPointWithVals(3, pointVal1);
    pointsArr[1] = createPointWithVals(3, pointVal2);
    pointsArr[2] = createPointWithVals(3, pointVal3);
    return pointsArr;
}

void Test1(bool isDebug) {
    Point **pointsArr;
    Matrix *W, *WA, *D1, *DA1, *D2, *DA2;

    /* Genarate points arr as input */
    pointsArr = pointsForTest1();
    W = computeMatrixW(pointsArr, 3);

    /* Calculate Matrix W */
    WA = createMatrix(3, 3, true);
    setMatrixValue(WA, 0 ,0, 0.0);
    setMatrixValue(WA, 0 ,1, exp(-1.5));
    setMatrixValue(WA, 0 ,2, exp(-6));
    setMatrixValue(WA, 1 ,1, 0.0);
    setMatrixValue(WA, 1 ,2, exp(-1.5));
    setMatrixValue(WA, 2 ,2, 0.0);

    if (isDebug == 1) {
        printf("points array: \n");
        printPointsArr(pointsArr, 3);
        printf("Matrix W calculated: \n");
        printMatrix(W);
        printf("Matrix A correct Matrix\n");
        printMatrix(WA);
    }
    
    (isMatrixEqual(W,WA)) ?
        printf("'test Calc Matrix W'\t\tresult: Great!\n") : 
        printf("'test Calc Matrix W'\t\tresult: Problem!\n");
    
    /* Calculate Matrix D */
    D1 = computeMatrixD(W);    
    DA1 = createMatrix(3, 3, true);
    setMatrixValue(DA1, 0 ,0, (0.0 + exp(-1.5) + exp(-6)));
    setMatrixValue(DA1, 1 ,1, (exp(-1.5) + 0 + exp(-1.5)));
    setMatrixValue(DA1, 2 ,2, (exp(-6) + exp(-1.5) + 0.0));

    if (isDebug == 1) {
        printf("Matrix W calculated: \n");
        printMatrix(W);
        printf("Matrix D calc:\n");
        printMatrix(D1);
        printf("Matrix A correct Matrix\n");
        printMatrix(DA1);
    }
    
    (isMatrixEqual(D1,DA1)) ?
        printf("'test Calc Matrix D'\t\tresult: Great!\n") : 
        printf("'test Calc Matrix D'\t\tresult: Problem!\n");

    /* Calculate Matrix D^-0.5 if Matrix D is good */
    if (isMatrixEqual(D1,DA1)) {
        D2 = computeMatrixDMinusHalf(D1);

        DA2 = createMatrix(3, 3, true);
        setMatrixValue(DA2, 0 ,0, 1 / sqrt(0.0 + exp(-1.5) + exp(-6)));
        setMatrixValue(DA2, 1 ,1, 1 / sqrt(exp(-1.5) + 0 + exp(-1.5)));
        setMatrixValue(DA2, 2 ,2, 1 / sqrt(exp(-6) + exp(-1.5) + 0.0));

        if (isDebug) {
            printf("Matrix D^-0.5 calculated: \n");
            printMatrix(D2);
            printf("Matrix A2 correct Matrix:\n");
            printMatrix(DA2);
        }
        
        (isMatrixEqual(D2,DA2)) ?
        printf("'test Calc Matrix D^-0.5'\tresult: Great!\n") : 
        printf("'test Calc Matrix D^-0.5'\tresult: Problem!\n");
    }

    freeMatrix(W);
    freeMatrix(WA);
    freeMatrix(D1);
    freeMatrix(D2);
    freeMatrix(DA1);
    freeMatrix(DA2);
    freeMemPointsArr(pointsArr, 3);
}

/*
void testCalcMatrixLnorm(bool isDebug) {
    check
}
*/

int main() {
    if (TestMode) {
        testMain(false);
    }
    return 0;
}