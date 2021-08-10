#include "spkmeans.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h> /* check */

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

Matrix* createUnitMatrix(int dim, bool isSymmetric) {
    Matrix* I; 
    int i;
    I = createMatrix(dim, dim, isSymmetric);
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
        data[i] = (double*) calloc(i + 1, sizeof(double));
        assert(data[i] != NULL);
    }

    return data;
}

Matrix* cloneMatrix(Matrix* A) {
    Matrix *Atag;
    int i, j;
    Atag = createMatrix(A -> rows, A -> cols, A -> isSymmetric);
    MatrixIterRows(A, i) {
        if (A -> isSymmetric){
            MatrixIterColsSym(A, i, j) {
                setMatrixValue(Atag, i, j, getMatrixValue(A,i,j));
            }
        }
        else {
            MatrixIterCols(A, j) {
                setMatrixValue(Atag, i, j, getMatrixValue(A,i,j));
            }
        }
    }
    return Atag;
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
                if (fabs(getMatrixValue(A,i,j) - getMatrixValue(B,i,j)) > epsilon) {
                    return false;
                }
            }
        }
        else {
            MatrixIterCols(A, j) {
                if (fabs(getMatrixValue(A,i,j) - getMatrixValue(B,i,j)) > epsilon) {
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

Point* createPointFromMatrixCol(Matrix* A, int col) {
    Point *p;
    int i;
    p = createPoint(A -> rows);
    MatrixIterRows(A, i) {
        setDataPointVal(p, i, getMatrixValue(A, i, col));
    }
    return p;
}

int compareEigens(const void * a, const void * b) {
    Eigen *A, *B;
    A = (Eigen *) a;
    B = (Eigen *) b;
    return (A -> value < B -> value) ? -1 : (A -> value > B -> value);
}

Eigen* getSortedEigen(Matrix* A) {
    Matrix *V;
    Eigen* eigens;
    int i;

    V = jacobiAlgo(&A);
    eigens = (Eigen*) calloc(V -> rows, sizeof(Eigen));
    assert(eigens != NULL);

    MatrixIterCols(V, i) {
        eigens[i].value = getMatrixValue(A, i, i);
        eigens[i].vector = createPointFromMatrixCol(V, i);
    }

    qsort(eigens, V -> rows, sizeof(Eigen), compareEigens);
    return eigens;
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
    double epsilon = 0.00000001;
    for (i = 0; i < point1 -> d; i++) {
        if (fabs(getDataPointVal(point1, i) - getDataPointVal(point2, i)) > epsilon) {
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
    return exp(-0.5 * dist);
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
    I = createUnitMatrix(W -> rows, true);
    tmp = multiply(multiply(D2, W), D2);
    return sub(I, tmp);
}

/* ################ */
/* Jacobi algorithm */
/* ################ */

MaxAbsulteValue getmaxAbsulteValue(Matrix* A) {
    int i, j;
    MaxAbsulteValue m;
    assert(A -> isSymmetric == true);

    m.value = 0;
    MatrixIterRows(A, i) {
        MatrixIterColsSym(A, i, j) {
            if(i == j) {
                continue;
            }
            if (fabs(getMatrixValue(A, i, j)) > fabs(m.value)) {
                m.i = j;
                m.j = i;
                m.value = getMatrixValue(A, i, j);
            }
        }
    }

    return m;
}

double getTheta(Matrix* A, MaxAbsulteValue mav) {
    double Aii, Aij, Ajj;

    Aij = mav.value;
    Aii = getMatrixValue(A, mav.i, mav.i);
    Ajj = getMatrixValue(A, mav.j, mav.j);

    return (Ajj - Aii) / (2 * Aij);
}

double getT(double theta) {
    int signTheta = theta >= 0 ? 1 : -1;
    return signTheta / (fabs(theta) + sqrt(theta * theta + 1));
}

double getC(double t) {
    return 1 / (sqrt(t * t + 1));
}

Matrix* createP(int dim, double c, double s, MaxAbsulteValue mav) {
    Matrix *P = createUnitMatrix(dim, false);
    setMatrixValue(P, mav.i, mav.i, c);
    setMatrixValue(P, mav.j, mav.j, c);
    if (mav.i > mav.j) {
        setMatrixValue(P, mav.i, mav.j, -s);
        setMatrixValue(P, mav.j, mav.i, s);
    } else {
        setMatrixValue(P, mav.i, mav.j, s);
        setMatrixValue(P, mav.j, mav.i, -s);
    }
    return P;
}

Matrix* createAtag(Matrix* A, double c, double s, MaxAbsulteValue mav) {
    int r;
    Matrix *Atag = cloneMatrix(A);

    MatrixIterRows(A, r) {
        if(r == mav.i || r == mav.j) {
            continue;
        }
        setMatrixValue(Atag, r, mav.i, c * getMatrixValue(A, r, mav.i) - s * getMatrixValue(A, r, mav.j));
        setMatrixValue(Atag, r, mav.j, c * getMatrixValue(A, r, mav.j) + s * getMatrixValue(A, r, mav.i));
    }
    setMatrixValue(Atag, mav.i, mav.i, c * c * getMatrixValue(A, mav.i, mav.i) + s * s * getMatrixValue(A, mav.j, mav.j) - 2 * s * c * getMatrixValue(A, mav.i, mav.j));
    setMatrixValue(Atag, mav.j, mav.j, s * s * getMatrixValue(A, mav.i, mav.i) + c * c * getMatrixValue(A, mav.j, mav.j) + 2 * s * c * getMatrixValue(A, mav.i, mav.j));
    setMatrixValue(Atag, mav.i, mav.j, 0);
    return Atag;
}

double getOffDiagMatrixSquareSum(Matrix* A) {
    int i, j;
    double sum = 0;
    assert(A -> isSymmetric);

    MatrixIterRows(A, i) {
        MatrixIterColsSym(A, i, j) {
            if(i == j) {
                continue;
            }
            sum += 2 * pow(getMatrixValue(A,i,j), 2);
        }
    }
    return sum;
}

bool isNeedToStopJabobi(Matrix* A, Matrix* Atag) {
    double epsilon = 0.001;
    return (fabs(getOffDiagMatrixSquareSum(A) - getOffDiagMatrixSquareSum(Atag)) <= epsilon);
}

Matrix* jacobiAlgo(Matrix** A_origin) {
    Matrix *Pi, *V, *PTemp, *Atag, *A;
    MaxAbsulteValue mav;
    bool isNeedToStop;
    double theta, c, s, t;
    A = *A_origin;
    assert(A -> rows == A -> cols);
    V = createUnitMatrix(A -> rows, false);

    while (true) {
        mav = getmaxAbsulteValue(A);
        if (mav.value == 0) { /*///////////////////////////////////////////////////////////need to talk about it///////////////////////////////////////////////////*/
            break;
        }
        theta = getTheta(A, mav);
        t = getT(theta);
        c = getC(t);
        s = t * c;
        Pi = createP(A -> rows, c, s, mav);
        PTemp = multiply(V, Pi);
        freeMatrix(V);
        V = PTemp;
        Atag = createAtag(A, c, s, mav);
        isNeedToStop = isNeedToStopJabobi(A, Atag);
        freeMatrix(A); /*//////////////////////////////????? maybe need original A in anther step? ////////////////////////////////////////////////////////////*/
        A = Atag;
        if (isNeedToStop) {
            break;
        }
    }
    *A_origin = A;
    return V;
}

/* ############# */
/* Tests section */
/* ############# */

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
    int dim = 3;
    double pointVal1[3] = {0,0,0};
    double pointVal2[3] = {1,1,1}; 
    double pointVal3[3] = {2,2,2};

    pointsArr = calloc(numOfPoints, sizeof(Point));
    assert(pointsArr != NULL); 
    
    pointsArr[0] = createPointWithVals(dim, pointVal1);
    pointsArr[1] = createPointWithVals(dim, pointVal2);
    pointsArr[2] = createPointWithVals(dim, pointVal3);
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
        printf("Test 1 - points array: \n");
        printPointsArr(pointsArr, 3);
        printf("Test 1 - Matrix W calculated: \n");
        printMatrix(W);
        printf("Test 1 - Matrix A correct Matrix\n");
        printMatrix(WA);
    }
    
    (isMatrixEqual(W,WA)) ?
        printf("Test1 - Matrix W\t\tresult: Great!\n") : 
        printf("Test1 - Matrix W\t\tresult: Problem!\n");
    
    /* Calculate Matrix D */
    D1 = computeMatrixD(W);    
    DA1 = createMatrix(3, 3, true);
    setMatrixValue(DA1, 0 ,0, (0.0 + exp(-1.5) + exp(-6)));
    setMatrixValue(DA1, 1 ,1, (exp(-1.5) + 0 + exp(-1.5)));
    setMatrixValue(DA1, 2 ,2, (exp(-6) + exp(-1.5) + 0.0));

    if (isDebug == 1) {
        printf("Test1 - Matrix W calculated: \n");
        printMatrix(W);
        printf("Test1 - Matrix D calc:\n");
        printMatrix(D1);
        printf("Test1 - Matrix A correct Matrix\n");
        printMatrix(DA1);
    }
    
    (isMatrixEqual(D1,DA1)) ?
        printf("Test1 - Matrix D\t\tresult: Great!\n") : 
        printf("Test1 - Matrix D\t\tresult: Problem!\n");

    /* Calculate Matrix D^-0.5 if Matrix D is good */
    if (isMatrixEqual(D1,DA1)) {
        D2 = computeMatrixDMinusHalf(D1);

        DA2 = createMatrix(3, 3, true);
        setMatrixValue(DA2, 0 ,0, 1 / sqrt(0.0 + exp(-1.5) + exp(-6)));
        setMatrixValue(DA2, 1 ,1, 1 / sqrt(exp(-1.5) + 0 + exp(-1.5)));
        setMatrixValue(DA2, 2 ,2, 1 / sqrt(exp(-6) + exp(-1.5) + 0.0));

        if (isDebug) {
            printf("Test1 - Matrix D^-0.5 calculated: \n");
            printMatrix(D2);
            printf("Test1 - Matrix A2 correct Matrix:\n");
            printMatrix(DA2);
        }
        
        (isMatrixEqual(D2,DA2)) ?
        printf("Test1 - Matrix D^-0.5\t\tresult: Great!\n") : 
        printf("Test1 - Matrix D^-0.5\t\tresult: Problem!\n");
    }

    freeMatrix(W);
    freeMatrix(WA);
    freeMatrix(D1);
    freeMatrix(D2);
    freeMatrix(DA1);
    freeMatrix(DA2);
    freeMemPointsArr(pointsArr, 3);
}

Point** pointsForTestE0() {
    Point **pointsArr;
    int numOfPoints = 5;
    int dim = 4;
    double pointVal1[4] = {0.1255,-0.4507,-0.232,-0.0987};
    double pointVal2[4] = {0.344,0.344,0.4419,-0.3662}; 
    double pointVal3[4] = {-0.1011,-0.2081,0.4794,-0.4699};
    double pointVal4[4] = {-0.3324,0.2877,0.3182,0.3166}; 
    double pointVal5[4] = {0.1958,-0.0248,0.0681,0.2088};

    pointsArr = calloc(numOfPoints, sizeof(Point));
    assert(pointsArr != NULL); 
    
    pointsArr[0] = createPointWithVals(dim, pointVal1);
    pointsArr[1] = createPointWithVals(dim, pointVal2);
    pointsArr[2] = createPointWithVals(dim, pointVal3);
    pointsArr[3] = createPointWithVals(dim, pointVal4);
    pointsArr[4] = createPointWithVals(dim, pointVal5);
    return pointsArr;
}

void TestE0(bool isDebug) {
    Point **pointsArr;
    Matrix *W, *WA, *D, *DA, *L, *LA;
    int numOfPoints = 5; 

    /* Genarate points arr as input */
    pointsArr = pointsForTestE0();
    W = computeMatrixW(pointsArr, numOfPoints);

    /* Calculate Matrix W */
    WA = createMatrix(numOfPoints, numOfPoints, true);
    setMatrixValue(WA, 0 ,0, 0.0);
    setMatrixValue(WA, 0 ,1, 0.5776);
    setMatrixValue(WA, 0 ,2, 0.6478);
    setMatrixValue(WA, 0 ,3, 0.5743);
    setMatrixValue(WA, 0 ,4, 0.7375);
    setMatrixValue(WA, 1 ,1, 0.0);
    setMatrixValue(WA, 1 ,2, 0.6985);
    setMatrixValue(WA, 1 ,3, 0.6155);
    setMatrixValue(WA, 1 ,4, 0.6728);
    setMatrixValue(WA, 2 ,2, 0.0);
    setMatrixValue(WA, 2 ,3, 0.6152);
    setMatrixValue(WA, 2 ,4, 0.6483);
    setMatrixValue(WA, 3 ,3, 0.0);
    setMatrixValue(WA, 3 ,4, 0.7148);
    setMatrixValue(WA, 4 ,4, 0.0);

    if (isDebug == 1) {
        printf("\nTestE0 - points array: \n");
        printPointsArr(pointsArr, numOfPoints);
        printf("\nTestE0 - Matrix W calculated: \n");
        printMatrix(W);
        printf("\nTestE0 - Matrix WA correct Matrix\n");
        printMatrix(WA);
    }
    
    (isMatrixEqual(W,WA)) ?
        printf("TestE0 - Matrix W\t\tresult: Great!\n") : 
        printf("TestE0 - Matrix W\t\tresult: Problem!\n");
    
    /* Calculate Matrix D */
    D = computeMatrixD(W);    
    
    DA = createMatrix(numOfPoints, numOfPoints, true);
    setMatrixValue(DA, 0 ,0, 2.5372);
    setMatrixValue(DA, 1 ,1, 2.5644);
    setMatrixValue(DA, 2 ,2, 2.6098);
    setMatrixValue(DA, 3 ,3, 2.5199);
    setMatrixValue(DA, 4 ,4, 2.7733);

    if (isDebug == 1) {
        printf("\nTestE0 - Matrix W calculated: \n");
        printMatrix(W);
        printf("\nTestE0 - Matrix D calc:\n");
        printMatrix(D);
        printf("\nTestE0 - Matrix A correct Matrix\n");
        printMatrix(DA);
    }
    
    (isMatrixEqual(D,DA)) ?
        printf("TestE0 - Matrix D\t\tresult: Great!\n") : 
        printf("TestE0 - Matrix D\t\tresult: Problem!\n");

    L = computeMatrixLnorm(W,D);
    
    LA = createMatrix(numOfPoints, numOfPoints, true);
    setMatrixValue(LA, 0 ,0, 1.0);
    setMatrixValue(LA, 0 ,1, -0.2264);
    setMatrixValue(LA, 0 ,2, -0.2517);
    setMatrixValue(LA, 0 ,3, -0.2271);
    setMatrixValue(LA, 0 ,4, -0.2780);

    setMatrixValue(LA, 1 ,1, 1.0);
    setMatrixValue(LA, 1 ,2, -0.2700);
    setMatrixValue(LA, 1 ,3, -0.2421);
    setMatrixValue(LA, 1 ,4, -0.2523);

    setMatrixValue(LA, 2 ,2, 1.0);
    setMatrixValue(LA, 2 ,3, -0.2399);
    setMatrixValue(LA, 2 ,4, -0.2410);

    setMatrixValue(LA, 3 ,3, 1.0);
    setMatrixValue(LA, 3 ,4, -0.2704);

    setMatrixValue(LA, 4 ,4, 1.0);

    if (isDebug == 1) {
        printf("\nTestE0 - Matrix L calculated: \n");
        printMatrix(L);
        printf("\nTestE0 - Matrix LA correct Matrix\n");
        printMatrix(LA);
    }
    
    (isMatrixEqual(L,LA)) ?
        printf("TestE0 - Matrix L\t\tresult: Great!\n") : 
        printf("TestE0 - Matrix L\t\tresult: Problem!\n");

    freeMatrix(W);
    freeMatrix(WA);
    freeMatrix(D);
    freeMatrix(DA);
    freeMatrix(L);
    freeMatrix(LA);
    freeMemPointsArr(pointsArr, 3);
}

Point** pointsForTestE1() {
    Point **pointsArr;
    int numOfPoints = 8;
    int dim = 5;
    double pointVal1[5] = {-0.1119,0.3605,0.2079,0.1336,0.0439};
    double pointVal2[5] = {-0.2852,0.3003,-0.0142,-0.0924,0.4535}; 
    double pointVal3[5] = {-0.1075,0.3295,0.4349,-0.4489,-0.4656};
    double pointVal4[5] = {-0.3084,0.1954,0.4023,-0.1842,0.0598}; 
    double pointVal5[5] = {0.378,0.0048,0.4656,-0.4093,0.2412};
    double pointVal6[5] = {-0.1625,0.1883,-0.0201,-0.0467,0.3151};
    double pointVal7[5] = {-0.4696,-0.2751,-0.4395,-0.3948,-0.0979};
    double pointVal8[5] = {-0.4219,0.4115,0.304,0.4548,0.1747};

    pointsArr = calloc(numOfPoints, sizeof(Point));
    assert(pointsArr != NULL); 
    
    pointsArr[0] = createPointWithVals(dim, pointVal1);
    pointsArr[1] = createPointWithVals(dim, pointVal2);
    pointsArr[2] = createPointWithVals(dim, pointVal3);
    pointsArr[3] = createPointWithVals(dim, pointVal4);
    pointsArr[4] = createPointWithVals(dim, pointVal5);
    pointsArr[5] = createPointWithVals(dim, pointVal6);
    pointsArr[6] = createPointWithVals(dim, pointVal7);
    pointsArr[7] = createPointWithVals(dim, pointVal8);
    return pointsArr;
}

void TestE1(bool isDebug) {
    Point **pointsArr;
    Matrix *W, *WA, *D, *DA, *L, *LA;
    int numOfPoints = 8; 

    /* Genarate points arr as input */
    pointsArr = pointsForTestE1();
    W = computeMatrixW(pointsArr, numOfPoints);

    /* Calculate Matrix W */
    WA = createMatrix(numOfPoints, numOfPoints, true);
    setMatrixValue(WA, 0 ,0, 0.0);
    setMatrixValue(WA, 0 ,1, 0.7598);
    setMatrixValue(WA, 0 ,2, 0.6679);
    setMatrixValue(WA, 0 ,3, 0.7975);
    setMatrixValue(WA, 0 ,4, 0.6455);
    setMatrixValue(WA, 0 ,5, 0.8041);
    setMatrixValue(WA, 0 ,6, 0.5717);
    setMatrixValue(WA, 0 ,7, 0.7875);

    setMatrixValue(WA, 1 ,1, 0.0);
    setMatrixValue(WA, 1 ,2, 0.5775);
    setMatrixValue(WA, 1 ,3, 0.7444);
    setMatrixValue(WA, 1 ,4, 0.6218);
    setMatrixValue(WA, 1 ,5, 0.8953);
    setMatrixValue(WA, 1 ,6, 0.6156);
    setMatrixValue(WA, 1 ,7, 0.6999);

    setMatrixValue(WA, 2 ,2, 0.0);
    setMatrixValue(WA, 2 ,3, 0.7273);
    setMatrixValue(WA, 2 ,4, 0.6318);
    setMatrixValue(WA, 2 ,5, 0.6063);
    setMatrixValue(WA, 2 ,6, 0.5535);
    setMatrixValue(WA, 2 ,7, 0.5594);

    setMatrixValue(WA, 3 ,3, 0.0);
    setMatrixValue(WA, 3 ,4, 0.6800);
    setMatrixValue(WA, 3 ,5, 0.7661);
    setMatrixValue(WA, 3 ,6, 0.6027);
    setMatrixValue(WA, 3 ,7, 0.7045);

    setMatrixValue(WA, 4 ,4, 0.0);
    setMatrixValue(WA, 4 ,5, 0.6584);
    setMatrixValue(WA, 4 ,6, 0.5180);
    setMatrixValue(WA, 4 ,7, 0.5331);

    setMatrixValue(WA, 5 ,5, 0.0);
    setMatrixValue(WA, 5 ,6, 0.6436);
    setMatrixValue(WA, 5 ,7, 0.7038);

    setMatrixValue(WA, 6 ,6, 0.0);
    setMatrixValue(WA, 6 ,7, 0.5091);

    setMatrixValue(WA, 7 ,7, 0.0);

    if (isDebug == 1) {
        printf("\nTestE1 - points array: \n");
        printPointsArr(pointsArr, numOfPoints);
        printf("\nTestE1 - Matrix W calculated: \n");
        printMatrix(W);
        printf("\nTestE1 - Matrix A correct Matrix\n");
        printMatrix(WA);
    }
    
    (isMatrixEqual(W,WA)) ?
        printf("TestE1 - Matrix W\t\tresult: Great!\n") : 
        printf("TestE1 - Matrix W\t\tresult: Problem!\n");
    
    /* Calculate Matrix D */
    D = computeMatrixD(W);    
    
    DA = createMatrix(numOfPoints, numOfPoints, true);
    setMatrixValue(DA, 0 ,0, 5.0340);
    setMatrixValue(DA, 1 ,1, 4.9143);
    setMatrixValue(DA, 2 ,2, 4.3239);
    setMatrixValue(DA, 3 ,3, 5.0225);
    setMatrixValue(DA, 4 ,4, 4.2886);
    setMatrixValue(DA, 5 ,5, 5.0778);
    setMatrixValue(DA, 6 ,6, 4.0143);
    setMatrixValue(DA, 7 ,7, 4.4974);

    if (isDebug == 1) {
        printf("\nTestE1 - Matrix W calculated: \n");
        printMatrix(W);
        printf("\nTestE1 - Matrix D calc:\n");
        printMatrix(D);
        printf("\nTestE1 - Matrix A correct Matrix\n");
        printMatrix(DA);
    }
    
    (isMatrixEqual(D,DA)) ?
        printf("TestE1 - Matrix D\t\tresult: Great!\n") : 
        printf("TestE1 - Matrix D\t\tresult: Problem!\n");

    L = computeMatrixLnorm(W,D);
    
    LA = createMatrix(numOfPoints, numOfPoints, true);
    setMatrixValue(LA, 0 ,0, 1.0);
    setMatrixValue(LA, 0 ,1, -0.1528);
    setMatrixValue(LA, 0 ,2, -0.1432);
    setMatrixValue(LA, 0 ,3, -0.1586);
    setMatrixValue(LA, 0 ,4, -0.1389);
    setMatrixValue(LA, 0 ,5, -0.1590);
    setMatrixValue(LA, 0 ,6, -0.1272);
    setMatrixValue(LA, 0 ,7, -0.1655);

    setMatrixValue(LA, 1 ,1, 1.0);
    setMatrixValue(LA, 1 ,2, -0.1253);
    setMatrixValue(LA, 1 ,3, -0.1498);
    setMatrixValue(LA, 1 ,4, -0.1354);
    setMatrixValue(LA, 1 ,5, -0.1792);
    setMatrixValue(LA, 1 ,6, -0.1386);
    setMatrixValue(LA, 1 ,7, -0.1489);

    setMatrixValue(LA, 2 ,2, 1.0);
    setMatrixValue(LA, 2 ,3, -0.1561);
    setMatrixValue(LA, 2 ,4, -0.1467);
    setMatrixValue(LA, 2 ,5, -0.1294);
    setMatrixValue(LA, 2 ,6, -0.1329);
    setMatrixValue(LA, 2 ,7, -0.1269);

    setMatrixValue(LA, 3 ,3, 1.0);
    setMatrixValue(LA, 3 ,4, -0.1465);
    setMatrixValue(LA, 3 ,5, -0.1517);
    setMatrixValue(LA, 3 ,6, -0.1342);
    setMatrixValue(LA, 3 ,7, -0.1482);

    setMatrixValue(LA, 4 ,4, 1.0);
    setMatrixValue(LA, 4 ,5, -0.1411);
    setMatrixValue(LA, 4 ,6, -0.1248);
    setMatrixValue(LA, 4 ,7, -0.1214);

    setMatrixValue(LA, 5 ,5, 1.0);
    setMatrixValue(LA, 5 ,6, -0.1426);
    setMatrixValue(LA, 5 ,7, -0.1473);

    setMatrixValue(LA, 6 ,6, 1.0);
    setMatrixValue(LA, 6 ,7, -0.1198);

    setMatrixValue(LA, 7 ,7, 1.0);

    if (isDebug == 1) {
        printf("\nTestE1 - Matrix L calculated: \n");
        printMatrix(L);
        printf("\nTestE1 - Matrix LA correct Matrix\n");
        printMatrix(LA);
    }
    
    (isMatrixEqual(L,LA)) ?
        printf("TestE1 - Matrix L\t\tresult: Great!\n") : 
        printf("TestE1 - Matrix L\t\tresult: Problem!\n");

    freeMatrix(W);
    freeMatrix(WA);
    freeMatrix(D);
    freeMatrix(DA);
    freeMatrix(L);
    freeMatrix(LA);
    freeMemPointsArr(pointsArr, numOfPoints);
}


void testJacobi(bool isDebug) {
    Matrix *V, *expectedV, *A, *expectedA;
    bool testResult;
    A = createMatrix(3, 3, true);
    setMatrixValue(A, 0, 0, 3.0);
    setMatrixValue(A, 1, 0, 2.0);
    setMatrixValue(A, 1, 1, 0.0);
    setMatrixValue(A, 2, 0, 4.0);
    setMatrixValue(A, 2, 1, 2.0);
    setMatrixValue(A, 2, 2, 3.0);

    expectedV = createMatrix(3, 3, false);
    setMatrixValue(expectedV, 0, 0, 1 / sqrt(2));
    setMatrixValue(expectedV, 0, 1, - 1 / (3 * sqrt(2)));
    setMatrixValue(expectedV, 0, 2, 2.0 / 3.0);
    setMatrixValue(expectedV, 1, 0, 0);
    setMatrixValue(expectedV, 1, 1, 4 / (3 * sqrt(2)));
    setMatrixValue(expectedV, 1, 2, 1.0 / 3.0);
    setMatrixValue(expectedV, 2, 0, - 1 / sqrt(2));
    setMatrixValue(expectedV, 2, 1, - 1 / (3 * sqrt(2)));
    setMatrixValue(expectedV, 2, 2, 2.0 / 3.0);

    expectedA = createMatrix(3, 3, true);
    setMatrixValue(expectedA, 0, 0, -1);
    setMatrixValue(expectedA, 1, 1, -1);
    setMatrixValue(expectedA, 2, 2, 8);

    V = jacobiAlgo(&A);
    testResult = isMatrixEqual(A, expectedA) && isMatrixEqual(V, expectedV);
    (testResult) ?
        printf("'test Jacobi'\t\t\tresult: Great!\n") :
        printf("'test Jacobi'\t\t\tresult: Problem!\n");
    if (isDebug == 1) {
        printf("Matrix A\n");
        printMatrix(A);
        printf("Expected Matrix A\n");
        printMatrix(expectedA);
        printf("Matrix V\n");
        printMatrix(V);
        printf("Expected Matrix V\n");
        printMatrix(expectedV);
    }
}

void testEigen(bool isDebug) {
    Matrix *A;
    bool testResult;
    Eigen *eigens, *expectedEigens;
    int i;
    A = createMatrix(3, 3, true);
    setMatrixValue(A, 0, 0, 3.0);
    setMatrixValue(A, 1, 0, 2.0);
    setMatrixValue(A, 1, 1, 0.0);
    setMatrixValue(A, 2, 0, 4.0);
    setMatrixValue(A, 2, 1, 2.0);
    setMatrixValue(A, 2, 2, 3.0);

    expectedEigens = (Eigen *) calloc(3, sizeof(Eigen));
    assert(expectedEigens != NULL);
    expectedEigens[0].value = -1;
    expectedEigens[0].vector = createPoint(3);
    expectedEigens[1].value = -1;
    expectedEigens[1].vector = createPoint(3);
    expectedEigens[2].value = 8;
    expectedEigens[2].vector = createPoint(3);
    setDataPointVal(expectedEigens[0].vector, 0, 1 / sqrt(2));
    setDataPointVal(expectedEigens[0].vector, 1, 0);
    setDataPointVal(expectedEigens[0].vector, 2, - 1 / sqrt(2));
    setDataPointVal(expectedEigens[1].vector, 0, - 1 / (3 * sqrt(2)));
    setDataPointVal(expectedEigens[1].vector, 1, 4 / (3 * sqrt(2)));
    setDataPointVal(expectedEigens[1].vector, 2, - 1 / (3 * sqrt(2)));
    setDataPointVal(expectedEigens[2].vector, 0, 2.0 / 3.0);
    setDataPointVal(expectedEigens[2].vector, 1, 1.0 / 3.0);
    setDataPointVal(expectedEigens[2].vector, 2, 2.0 / 3.0);

    eigens = getSortedEigen(A);
    testResult = true;
    for(i = 0; i < 3; i++) {
        if (eigens[i].value != expectedEigens[i].value) {
            testResult = false;
            break;
        }
        if(!isPointsEquel(eigens[i].vector, expectedEigens[i].vector)) {
            testResult = false;
            break;
        }
    }
    (testResult) ?
        printf("'test Eigens'\t\t\tresult: Great!\n") :
        printf("'test Eigens'\t\t\tresult: Problem!\n");
    if (isDebug == 1) {
        printf("Eigens\n");
    }
}

void testMain(bool isDebug) {
    testMultiplyMatrixs(isDebug);
    /*Test1(isDebug);*/
    testJacobi(isDebug);
    testEigen(isDebug);
    TestE0(1);
    TestE1(isDebug);
}

int main() {
    int i;
    double arr1[4] = {0.1255,-0.4507,-0.2320,-0.0987}; 
    double arr2[4] = {0.3440,0.3440,0.4419,-0.3662};
    double tmp, res = 0;
    for (i = 0; i < 4; i++) {
        tmp = arr1[i] - arr2[i]; 
        res += tmp * tmp;
    }
    res = sqrt(res); 
    res = res*-0.5; 
    res = exp(res);
    printf("res = %f", res);

    if (TestMode) {
        testMain(false);
    }

    return 0;
}