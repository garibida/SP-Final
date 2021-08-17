#include "spkmeans.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* ################################################################################################ */
/*                               Matrix operations section                                          */
/* ################################################################################################ */

Matrix* createMatrix(int rows, int cols, bool isSymmetric) {
    Matrix* A = (Matrix*) malloc(sizeof(Matrix));
    assert(A != NULL);
    A -> rows = rows;
    A -> cols = cols;
    A -> isSymmetric = isSymmetric;
    A -> data = isSymmetric ? createSymmetricMatrixData(rows) : createMatrixData(rows, cols);
    return A;
}

 /* ################################################################################################ */

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

 /* ################################################################################################ */

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

 /* ################################################################################################ */

Matrix* createUnitMatrix(int dim, bool isSymmetric) {
    Matrix* I; 
    int i;
    I = createMatrix(dim, dim, isSymmetric);
    MatrixIterRows(I, i) {
        setMatrixValue(I, i, i, 1.0);
    }
    return I;
}

 /* ################################################################################################ */

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

 /* ################################################################################################ */

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

 /* ################################################################################################ */

double getMatrixValue(Matrix* A, int row, int col) {
    assert(row < (A -> rows) && col < (A -> cols));
    return ((A -> isSymmetric) && col > row) ?
                (A -> data)[col][row]:
                (A -> data)[row][col];
}

 /* ################################################################################################ */

void setMatrixValue(Matrix* A, int row, int col, double value) {
    assert(row < (A -> rows) && col < (A -> cols));
    if ((A -> isSymmetric) && col > row) {
        (A -> data)[col][row] = value;
    } else {
        (A -> data)[row][col] = value;
    }
}

 /* ################################################################################################ */

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

 /* ################################################################################################ */

Matrix* add(Matrix *A, Matrix *B) {
    Matrix* C;
    int i, j;
    bool isSymmetric;
    assert(A -> rows == B -> rows);
    assert(A -> cols == B -> cols);

    isSymmetric = (A -> isSymmetric) && (B -> isSymmetric);

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

 /* ################################################################################################ */

Matrix* sub(Matrix *A, Matrix *B) {
    Matrix* C;
    int i, j;
    bool isSymmetric;
    assert(A -> rows == B -> rows);
    assert(A -> cols == B -> cols);

    isSymmetric = (A -> isSymmetric) && (B -> isSymmetric);

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

 /* ################################################################################################ */

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

 /* ################################################################################################ */

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
    bool isSymmetric;
    assert(A -> rows == B -> rows);
    assert(A -> cols == B -> cols);

    isSymmetric = (A -> isSymmetric) && (B -> isSymmetric);

    MatrixIterRows(A, i) {
        if (isSymmetric) {
            MatrixIterColsSym(A, i, j) {
                if (fabs(getMatrixValue(A,i,j) - getMatrixValue(B,i,j)) > EPSILON) {
                    return false;
                }
            }
        }
        else {
            MatrixIterCols(A, j) {
                if (fabs(getMatrixValue(A,i,j) - getMatrixValue(B,i,j)) > EPSILON) {
                    return false;
                }
            }
        }
    }
    return true;
}

 /* ################################################################################################ */

void printMatrix(Matrix* A) { /* for debug ################################################################################## */
    int i, j;
    printf("===================\n");
    MatrixIterRows(A,i) {
        MatrixIterCols(A,j) {
            printf("%.4f", getMatrixValue(A,i,j));
            if(j != (A->cols) - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

 /* ################################################################################################ */

Point* createPointFromMatrixCol(Matrix* A, int col) {
    Point *point;
    int i;
    point = createPoint(A -> rows);
    MatrixIterRows(A, i) {
        setDataPointVal(point, i, getMatrixValue(A, i, col));
    }
    return point;
}

 /* ################################################################################################ */

Point* createPointFromMatrixRow(Matrix* A, int row) {
    Point *p;
    int i;
    p = createPoint(A -> cols);
    MatrixIterCols(A, i) {
        setDataPointVal(p, i, getMatrixValue(A, row, i));
    }
    return p;
}

 /* ################################################################################################ */

int compareEigens(const void *a, const void *b) {
    Eigen *A, *B;
    A = (Eigen *) a;
    B = (Eigen *) b;
    return (A -> value < B -> value) ? -1 : (A -> value > B -> value);
}

 /* ################################################################################################ */

Eigens_Arr* getSortedEigen(Matrix *A) {
    Matrix *V;
    Eigens_Arr *eigens;
    int i;

    eigens = (Eigens_Arr*) malloc(sizeof(Eigens_Arr));
    assert(eigens != NULL);

    V = jacobiAlgo(&A);
    eigens->length = V->rows;
    eigens->arr = (Eigen*) calloc(eigens->length, sizeof(Eigen));
    assert(eigens->arr != NULL);

    MatrixIterCols(V, i) {
        (eigens->arr)[i].value = getMatrixValue(A, i, i);
        (eigens->arr)[i].vector = createPointFromMatrixCol(V, i);
    }

    qsort(eigens->arr, eigens->length, sizeof(Eigen), compareEigens); 
    /* ########################################################### check if the in order of vector of the same value is meaningful */ 
    return eigens;
}

 /* ################################################################################################ */

void freeEigens(Eigens_Arr *eigens) {
    int i;
    for (i = 0; i < (eigens->length); i++) {
        freeMemPoint((eigens->arr)[i].vector);
    }
    free(eigens);    
}

 /* ################################################################################################ */

PointsArray* matrixToPointsArray(Matrix *A) {
    PointsArray *points;
    Point *point;
    int i;

    points = createPointsArr(A->rows);
    MatrixIterRows(A, i) {
        point = createPointFromMatrixRow(A, i);
        setPointInArr(points, i, point);
    }

    return points;
}

/* ################################################################################################ */
/*                               Points operation section                                           */
/* ################################################################################################ */

Point* createPoint(int d) {
    Point *point;
    Point_data data;
    point = (Point*) malloc(sizeof(Point));
    assert(point != NULL);
    point->d = d;
    data = (Point_data) calloc(d, sizeof(double));
    assert(data != NULL);
    point->data = data;
    return point;
}

 /* ################################################################################################ */

Point* createPointWithVals(int d, double *values) {
    Point *point;
    int i;    
    point = createPoint(d);
    for (i = 0; i < d; i++) {
        setDataPointVal(point, i, values[i]);
    }
    return point;
}

 /* ################################################################################################ */

void setDataPointVal(Point *point, int index, double value) {
    (point->data)[index] = value;
}

 /* ################################################################################################ */

double getDataPointVal(Point *point, int index) {
    return point->data[index];
}

 /* ################################################################################################ */

void printPoint(Point *point) {
    int i, dim;
    dim = point->d;
    for (i = 0; i < dim; i++) {
        printf("%.4f", point -> data[i]);
        if (i != dim - 1) {
            printf(",");
        }
    }
    printf("\n");
}

 /* ################################################################################################ */

int isPointsEquel(Point *point1, Point *point2){
    int i;
    for (i = 0; i < point1 -> d; i++) {
        if (fabs(getDataPointVal(point1, i) - getDataPointVal(point2, i)) > EPSILON) {
            return false;
        }
    }
    return true;
} 

 /* ################################################################################################ */

void freeMemPoint(Point *point) {
    if (point != NULL) {
        free(point -> data);
        free(point); 
    }
}

 /* ################################################################################################ */

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

/* ################################################################################################ */
/*                               Points Array operation section                                     */
/* ################################################################################################ */

PointsArray* createPointsArr(int n)  {
    PointsArray* pointsArr = (PointsArray*) malloc(sizeof(pointsArr));
    assert(pointsArr != NULL);
    pointsArr->n = n;
    pointsArr->points = (Point**) calloc(n, sizeof(Point *));
    return pointsArr;
}

 /* ################################################################################################ */

Point* getPointFromArr(PointsArray* pointsArr, int i) {
    assert(i < pointsArr->n);
    return (pointsArr->points)[i];
}

 /* ################################################################################################ */

void setPointInArr(PointsArray* pointsArr, int i, Point* point) {
    assert(i < pointsArr->n);
    (pointsArr->points)[i] = point;
}

 /* ################################################################################################ */

void reallocPointsArr(PointsArray* pointsArr, int n) {
    pointsArr->points = (Point **) realloc(pointsArr->points, n * sizeof(Point*));
    assert(pointsArr != NULL);
    pointsArr->n = n;
}

 /* ################################################################################################ */

void printPointsArr(PointsArray *pointsArr) { /* for debug ####################################################################### */
    int i;
    for (i = 0; i < (pointsArr->n); i++) {
        printPoint(getPointFromArr(pointsArr,i));
    }
}

 /* ################################################################################################ */

void freeMemPointsArr(PointsArray *pointsArr) {
    int i;
    for (i = 0; i < (pointsArr->n); i++) {
        freeMemPoint(getPointFromArr(pointsArr, i));
    }
    free(pointsArr->points);
    free(pointsArr);
}

/* ################################################################################################ */
/*                                       Algorithm                                                  */
/* ################################################################################################ */

double computeDist(Point *point1, Point *point2) {
    double dist = 0, tmp = 0; 
    int i;
    assert(point1->d == point2->d);
    for (i = 0; i < point1->d; i++) {
        tmp = getDataPointVal(point1, i) - getDataPointVal(point2,i);
        dist += tmp * tmp;
    }
    return sqrt(dist);
}

 /* ################################################################################################ */

double computeDistW(Point *point1, Point *point2) {
    double dist = computeDist(point1, point2);
    return exp(-0.5 * dist);
}

 /* ################################################################################################ */

Matrix* computeMatrixW(PointsArray *pointsArr) {
    Matrix *W; 
    Point *pointI, *pointJ;
    int i, j;
    double wVal;
    W = createMatrix(pointsArr->n, pointsArr->n, true);
    
    MatrixIterRows(W, i) {
        pointI = getPointFromArr(pointsArr, i);
        MatrixIterColsSym(W, i, j) {
            pointJ = getPointFromArr(pointsArr, j);
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

 /* ################################################################################################ */

Matrix* computeMatrixD(Matrix *W) {
    Matrix *D; 
    int i, j, dim;
    double dVal = 0;
    dim = W->cols;
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

 /* ################################################################################################ */

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

 /* ################################################################################################ */

Matrix* computeMatrixL(Matrix *W, Matrix *D) { /* delete??? #################################################################### */
    return sub(D, W);
}

 /* ################################################################################################ */

Matrix* computeMatrixLnorm(Matrix *W ,Matrix *D) {
    Matrix *D2, *I, *tmp;
    D2 = computeMatrixDMinusHalf(D);
    I = createUnitMatrix(W->rows, true);
    tmp = multiply(multiply(D2, W), D2);
    return sub(I, tmp);
}

 /* ################################################################################################ */

int eigengapGetK(Eigens_Arr* eigens) {
    int delta, max_delta = 0, max_index = 0, i;
    Eigen *arr = eigens->arr;

    for (i = 0; i < (eigens->length) / 2; i++) {
        assert(arr[i].value <= arr[i + 1].value); /* ########################################################### FOR DEBUG */
        delta = fabs(arr[i].value - arr[i + 1].value);
        if (delta > max_delta) {
            max_delta = delta;
            max_index = i;
        }
    }

    return max_index + 1;
}

 /* ################################################################################################ */

Matrix* computeMatrixU(Eigens_Arr* eigens, int k) {
    Matrix* U;
    int i, j;
    U = createMatrix((eigens->arr)[0].vector->d, k, false);
    MatrixIterRows(U, i) {
        MatrixIterCols(U, j) {
            setMatrixValue(U, i, j, ((eigens->arr)[j].vector->data)[i]);
        }
    }
    return U;
}

 /* ################################################################################################ */

double* getRowsSqureRootSum(Matrix* U) {
    int i, j;
    double *squreSumPerCol;

    squreSumPerCol = (double*) calloc(U->rows, sizeof(double));
    assert(squreSumPerCol != NULL);

    MatrixIterRows(U, i) {
        MatrixIterCols(U, j) {
            squreSumPerCol[i] += pow(getMatrixValue(U, i, j), 2);
        }
    }

    MatrixIterRows(U, i) {
        squreSumPerCol[i] = sqrt(squreSumPerCol[i]);
    }

    return squreSumPerCol;
}

 /* ################################################################################################ */

Matrix* computeMatrixT(Matrix* U) {
    Matrix *T;
    int i, j;
    double *squreSumPerRow = getRowsSqureRootSum(U);
    T = createMatrix(U->rows, U->cols, false);

    MatrixIterRows(U, i) {
        MatrixIterCols(U, j) {
            setMatrixValue(T, i, j, getMatrixValue(U, i, j) / squreSumPerRow[i]);
        }
    }

    free(squreSumPerRow);
    return T;
}

/* ################################################################################################ */
/*                                       Jacobi algorithm                                           */
/* ################################################################################################ */

MaxAbsulteValue getMaxAbsulteValue(Matrix* A) {
    int i, j;
    MaxAbsulteValue m;
    assert(A->isSymmetric == true);

    m.value = 0;
    MatrixIterRows(A, i) {
        MatrixIterColsSym(A, i, j) {
            if (i == j) {
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

 /* ################################################################################################ */

double getTheta(Matrix* A, MaxAbsulteValue mav) {
    double Aii, Aij, Ajj;
    Aij = mav.value;
    Aii = getMatrixValue(A, mav.i, mav.i);
    Ajj = getMatrixValue(A, mav.j, mav.j);
    return (Ajj - Aii) / (2 * Aij);
}

 /* ################################################################################################ */

double getT(double theta) {
    int signTheta = (theta >= 0) ? 1 : -1;
    return signTheta / (fabs(theta) + sqrt(theta * theta + 1));
}

 /* ################################################################################################ */

double getC(double t) {
    return 1 / (sqrt(t * t + 1));
}

 /* ################################################################################################ */

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

 /* ################################################################################################ */

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

 /* ################################################################################################ */

double getOffDiagMatrixSquareSum(Matrix* A) {
    int i, j;
    double sum = 0;
    assert(A->isSymmetric);

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

 /* ################################################################################################ */

bool isNeedToStopJacobi(Matrix* A, Matrix* Atag) {
    return (fabs(getOffDiagMatrixSquareSum(A) - getOffDiagMatrixSquareSum(Atag)) <= EPSILON_YACOBI); 
}

 /* ################################################################################################ */

void calcJacobiParams(Matrix* A, MaxAbsulteValue mav, double *c, double *s) {
    double theta, t;
    theta = getTheta(A, mav);
    t = getT(theta);
    *c = getC(t);
    *s = t * (*c);
}

 /* ################################################################################################ */

void calcJacobiV(Matrix* A, MaxAbsulteValue mav, double c, double s, Matrix** V) {
    Matrix *Pi, *PTemp;
    Pi = createP(A -> rows, c, s, mav);
    PTemp = multiply(*V, Pi);
    freeMatrix(*V);
    *V = PTemp;
}

 /* ################################################################################################ */

Matrix* jacobiAlgo(Matrix** A_origin) {
    Matrix *V, *Atag, *A;
    MaxAbsulteValue mav;
    bool isNeedToStop;
    int i;
    const int MAX_ITER = 100;
    double c, s;
    A = *A_origin;
    assert(A -> rows == A -> cols);
    V = createUnitMatrix(A -> rows, false);

    for (i = 0; i < MAX_ITER; i++) {
        mav = getMaxAbsulteValue(A);
        if (mav.value == 0) { /* the Matrix is diagonal */ 
            break;
        }

        calcJacobiParams(A, mav, &c, &s);
        calcJacobiV(A, mav, c, s, &V);
        Atag = createAtag(A, c, s, mav);
        isNeedToStop = isNeedToStopJacobi(A, Atag);
        freeMatrix(A); /* ################################################################################# maybe need original A in anther step? ////////////////////////////////////////////////////////////*/
        A = Atag;
        if (isNeedToStop) {
            break;
        }
    }
    *A_origin = A;
    return V;
}

/* ################################################################################################ */
/*                                       Link list operations                                       */
/* ################################################################################################ */

void addToList(linked_list* list, Point* point) {
    node *n = (node*)malloc(sizeof(node));
    assert(n != NULL);
    n -> point = point;
    n -> next = NULL;
    (list->length)++;
    if(list -> head == NULL) {
        list -> head = n;
        list -> tail = n;
    } else {
        list -> tail -> next = n;
        list -> tail = n;
    }
}

 /* ################################################################################################ */

void freeList(linked_list* list, int isDeletePoint) {
    freeNode(list -> head, isDeletePoint);
    free(list);
}

 /* ################################################################################################ */

void freeNode(node* n, int isDeletePoint) {
    if (n != NULL) {
        freeNode(n -> next, isDeletePoint);
        if(isDeletePoint == true){
            free(n -> point);
        }
        free(n);
    }
}


/* ################################################################################################ */
/*                                       K-means algorithm                                          */
/* ################################################################################################ */

PointsArray* kmeans(PointsArray *pointsArr, PointsArray *centroidsArr, int k, int max_iter) {
    int iter, isChanged;

    for (iter = 0; iter < max_iter; iter++) {
        isChanged = computeCluster(k, centroidsArr, pointsArr);
        if (!isChanged) {
            break;
        }
    }
    return centroidsArr;
}

 /* ################################################################################################ */

bool computeCluster(int k, PointsArray *centroidsArr, PointsArray *pointsArr) {
    double minDist, dist;
    int minIndex, i, j;
    bool isChanged;
    Point* point;
    linked_list** clusters = (linked_list**)calloc(k, sizeof(linked_list*));
    assert(clusters != NULL);

    for(i = 0; i < k; i++) {
        clusters[i] = (linked_list*)calloc(1, sizeof(linked_list));
        assert(clusters[i] != NULL);
    }

    for (i = 0; i < (pointsArr->n); i++) {
        point = getPointFromArr(pointsArr, i);
        minIndex = 0;
        minDist = computeDist(getPointFromArr(centroidsArr, 0), point);
        for (j = 1; j < k; j++) {
            dist = computeDist(getPointFromArr(centroidsArr, j), point);
            if (dist < minDist) {
                minDist = dist;
                minIndex = j;
            }
        }
        addToList(clusters[minIndex], point);
    }

    isChanged = computeNewCentroids(clusters, centroidsArr, k);
    for (i = 0; i < k; i++){
        freeList(clusters[i], false);
    }
    free(clusters);
    return isChanged;
}

 /* ################################################################################################ */

bool computeNewCentroids(linked_list** clusters, PointsArray *centroidsArr, int k) {
    Point* oldCentroid, *newCentroid;
    int i, j, t, isChanged = false;
    double temp;
    node* n;
    for (i = 0; i < k; i++) {
        oldCentroid = copy_point(getPointFromArr(centroidsArr, i));
        newCentroid = getPointFromArr(centroidsArr, i);
        j = 0;
        for (n = (clusters[i]) -> head; n != NULL; n = n-> next) {
            for (t = 0; t < newCentroid->d; t++) {
                temp = getDataPointVal(newCentroid, t);
                temp = (temp * j + getDataPointVal(n->point,t)) / (j + 1);
                setDataPointVal(newCentroid, t, temp);
            }
            j++;
        }
        if (isPointsEquel(oldCentroid, newCentroid) == false) {
            isChanged = true;
        }
        free(oldCentroid);
    }
    return isChanged;
}

 /* ################################################################################################ */

PointsArray* getIntialCentroids(PointsArray *pointsArr, int k) {
    PointsArray *centroidsArr;
    int i;

    centroidsArr = createPointsArr(k);
    for (i = 0; i < k; i++) {
        setPointInArr(centroidsArr, i, copy_point(getPointFromArr(pointsArr, i)));
    }

    return centroidsArr;
}

 /* ################################################################################################ */

void printCentroids(PointsArray* centroids) {
    int i, j;
    Point *centroid;
    for (i = 0; i < (centroids->n); i++) {
        centroid = getPointFromArr(centroids, i);
        for (j = 0; j < (centroid->d); j++) {
            printf("%.4f", (centroid->data)[j]);
            if(j != (centroid->d) - 1){
                printf(",");
            }
        }
        printf("\n");
    }
}

/* ################################################################################################ */
/*                                       Main Functions                                             */
/* ################################################################################################ */

void matrixPrinter(PointsArray *points ,Goal goal) {
    Matrix *W, *D, *Lnorm;
    Eigens_Arr *eigens;

    W = computeMatrixW(points);
    freeMemPointsArr(points); /* #############################################################################need later?########################################################### */
    if (goal == wam) {
        printMatrix(W); /* #################################################################Check printing##################################################### */
        freeMatrix(W);
        return;
    }

    D = computeMatrixD(W);
    if (goal == ddg) {
        printMatrix(D); /* #################################################################Check printing##################################################### */
        freeMatrix(W);
        freeMatrix(D);
        return;
    }

    Lnorm = computeMatrixLnorm(W, D);
    freeMatrix(W);
    freeMatrix(D);
    if (goal == lnorm) {
        printMatrix(Lnorm); /* #################################################################Check printing##################################################### */
        freeMatrix(Lnorm);
        return;
    }

    eigens = getSortedEigen(Lnorm);
    freeMatrix(Lnorm);
    if (goal == jacobi) {
        /* #############################################################################print Eigens########################################################### */
    }
    freeEigens(eigens);
}

 /* ################################################################################################ */

int doSpk(PointsArray **points, int k) {
    Matrix *W, *D, *Lnorm, *U, *T;
    Eigens_Arr *eigens;

    /* Stage 1 */
    W = computeMatrixW(*points);
    freeMemPointsArr(*points); /* #############################################################################need later?########################################################### */
    /* Stage 2 */
    D = computeMatrixD(W);
    Lnorm = computeMatrixLnorm(W, D);
    freeMatrix(W);
    freeMatrix(D);
    /* Stage 3 */
    eigens = getSortedEigen(Lnorm);
    freeMatrix(Lnorm);
    k = (k == 0) ? eigengapGetK(eigens) : k;
    /* Stage 4 */
    U = computeMatrixU(eigens, k);
    freeEigens(eigens);
    /* Stage 5 */
    T = computeMatrixT(U);
    freeMatrix(U);

    *points = matrixToPointsArray(T);
    freeMatrix(T);
    return k;
}

/* ################################################################################################ */
/*                                       Get input and validation                                   */
/* ################################################################################################ */

PointsArray* readPointsArray(char *path) {
    double value, *firstPointValues;
    int d = 0, numOfPoints = 0;
    char ch;
    Point *point;
    PointsArray *pointsArr;
    FILE *input;
    int i = 0;

    pointsArr = createPointsArr(MAX_NUMBER_OF_POINTS); /* free mem */

    firstPointValues = (double*)calloc(MAX_FEATURES, sizeof(double));
    assert(firstPointValues != NULL);

    /* read points from file */
    input = fopen(path, "r");
    assert(input != NULL);

    while ((!feof(input)) && (d < MAX_FEATURES) ) {
        fscanf(input, "%lf%c", &value, &ch);
        firstPointValues[d] = value;
        d++;
        if (ch == '\n' || ch == '\r') {
            break;
        }
    }

    setPointInArr(pointsArr, 0, createPointWithVals(d, firstPointValues));
    numOfPoints++;
    free(firstPointValues);
    point = createPoint(d);

    while(!feof(input)) {
        fscanf(input, "%lf%c", &value, &ch);
        setDataPointVal(point, i ,value);
        i++;
        if (i == d) {
            setPointInArr(pointsArr, numOfPoints, point);
            point = createPoint(d);
            numOfPoints++;
            i = 0;
        }
    }
    
    if (numOfPoints < MAX_NUMBER_OF_POINTS) {
        reallocPointsArr(pointsArr, numOfPoints);
    }
    freeMemPoint(point);
    fclose(input);
    return pointsArr;
}

 /* ################################################################################################ */

Goal decide_command(char *arg) {
    int enumIndex;
    char *commands[] = {"spk", "wam", "ddg", "lnorm", "jacobi"};
    for (enumIndex = 0; enumIndex < ENUM_COUNT; enumIndex++) {
        if (strcmp(commands[enumIndex], arg) == 0) {
            return enumIndex; /* since spk == 0 in enum Goal, the correct value will set, if found */
        }
    }
    printf("Invalid Input!\n"); /* "%s is not a goal.\nchoose from: spk / wam / ddg / lnorm / yacobi\nexits...\n", arg */ 
    assert(0);
    return 0;
}

 /* ################################################################################################ */

int main(int argc, char *argv[]) {
    PointsArray *points, *centroids;
    Goal goal;
    char *path;
    int k, max_iter = 100; /*################################################WHAT???################################## */
    assert( !(argc == 3) ); /* if k not provided set to 0 or exit? */ 
    
    path = argv[3];
    points = readPointsArray(path);
    
    k = atoi(argv[1]);
    if ((k < 0) || (k > points->n)) {
        printf("Invalid Input!\n");
        assert(0);
    } 

    goal = decide_command(argv[2]);

    if (goal != spk) {
        matrixPrinter(points, goal);
        return 0; /* check how to exit */ 
    }
    
    k = doSpk(&points, k);
    centroids = getIntialCentroids(points, k);
    centroids = kmeans(points, centroids, k, max_iter);
    printCentroids(centroids);
    freeMemPointsArr(centroids);
    return 0; /* check how to exit */ 
}