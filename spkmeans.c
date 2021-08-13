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

Matrix* createUnitMatrix(int dim, bool isSymmetric) {
    Matrix* I; 
    int i;
    I = createMatrix(dim, dim, isSymmetric);
    MatrixIterRows(I, i) {
        setMatrixValue(I, i, i, 1.0);
    }
    return I;
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
    assert(row < (A -> rows) && col < (A -> cols));
    return ((A -> isSymmetric) && col > row) ?
                (A -> data)[col][row]:
                (A -> data)[row][col];
}

void setMatrixValue(Matrix* A, int row, int col, double value) {
    assert(row < (A -> rows) && col < (A -> cols));
    if ((A -> isSymmetric) && col > row) {
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
    updateMatrixSymmertircStatus(C); /* ########################################################### ask Garibi if needed */
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
    Point *point;
    int i;
    point = createPoint(A -> rows);
    MatrixIterRows(A, i) {
        setDataPointVal(point, i, getMatrixValue(A, i, col));
    }
    return point;
}

int compareEigens(const void *a, const void *b) {
    Eigen *A, *B;
    A = (Eigen *) a;
    B = (Eigen *) b;
    return (A -> value < B -> value) ? -1 : (A -> value > B -> value);
}

Eigens_Arr* getSortedEigen(Matrix* A) {
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

    qsort(eigens->arr, eigens->length, sizeof(Eigen), compareEigens); /* ########################################################### check if the in order of vector of the same value is meaningful */ 
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
    point->d = d;
    data = (Point_data) calloc(d, sizeof(double));
    assert(data != NULL);
    point->data = data;
    return point;
}

Point* createPointWithVals(int d, double *values) {
    Point *point;
    int i;    
    point = createPoint(d);
    for (i = 0; i < d; i++) {
        setDataPointVal(point, i, values[i]);
    }
    return point;
}

void setDataPointVal(Point *point, int index, double value) {
    (point->data)[index] = value;
}

double getDataPointVal(Point *point, int index) {
    return point->data[index];
}

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

void printPointsArr(Point **pointArr, int numOfPoints) {
    int i;
    for (i = 0; i < numOfPoints; i++) {
        printPoint(pointArr[i]);
    }
}

int isPointsEquel(Point *point1, Point *point2){
    int i;
    for (i = 0; i < point1 -> d; i++) {
        if (fabs(getDataPointVal(point1, i) - getDataPointVal(point2, i)) > EPSILON) {
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

int eigengapGetK(Eigens_Arr* eigens) {
    int delta, max_delta = 0, max_index = 0, i;
    Eigen *arr = eigens -> arr;

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

/* ################ */
/* Jacobi algorithm */
/* ################ */

MaxAbsulteValue getMaxAbsulteValue(Matrix* A) {
    int i, j;
    MaxAbsulteValue m;
    assert(A -> isSymmetric == true);

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

double getTheta(Matrix* A, MaxAbsulteValue mav) {
    double Aii, Aij, Ajj;
    Aij = mav.value;
    Aii = getMatrixValue(A, mav.i, mav.i);
    Ajj = getMatrixValue(A, mav.j, mav.j);
    return (Ajj - Aii) / (2 * Aij);
}

double getT(double theta) {
    int signTheta = (theta >= 0) ? 1 : -1;
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
    return (fabs(getOffDiagMatrixSquareSum(A) - getOffDiagMatrixSquareSum(Atag)) <= EPSILON); 
    /* ################################################################################## ask Garibi - is 0.0001 ok or should it stay 0.001? */ 
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
        mav = getMaxAbsulteValue(A);
        if (mav.value == 0) { /* the Matrix is diagonal */ 
            break;
        }

        /* ################################################################################# ask Garibi - maybe we should calc all params in a function for readability? */ 
        /* something like calcParams(&theth, &t, &c, &s) */ 

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
        freeMatrix(A); /* ################################################################################# maybe need original A in anther step? ////////////////////////////////////////////////////////////*/
        A = Atag;
        if (isNeedToStop) {
            break;
        }
    }
    *A_origin = A;
    return V;
}

/* ######################## */
/* Get input and validation */
/* ######################## */

Point** readPointsArray(char *path, int *d, int *numOfPoints) {
    double value, *firstPointValues;
    char ch;
    Point **pointsArr, *point;
    FILE *input;
    int i = 0;
    *d = 0; 
    *numOfPoints = 0;

    pointsArr = (Point**)calloc(MAX_NUMBER_OF_POINTS, sizeof(Point*)); /* free mem */ 
    assert(pointsArr != NULL);

    firstPointValues = (double*)calloc(MAX_FEATURES, sizeof(double));
    assert(firstPointValues != NULL);

    /* read points from file */
    input = fopen(path, "r");
    assert(input != NULL);

    while ( ( !feof(input) ) && ( (*d) < MAX_FEATURES) ) {
        fscanf(input, "%lf%c", &value, &ch);
        firstPointValues[ (*d)++ ] = value;
        if (ch == '\n') {
            break;
        }
    }

    pointsArr[0] = createPointWithVals( (*d), firstPointValues);
    (*numOfPoints)++;
    free(firstPointValues);
    point = createPoint( (*d) );  

    while(!feof(input)) {
        fscanf(input, "%lf%c", &value, &ch);
        setDataPointVal(point, i ,value);
        i++;
        if (i == (*d) ) { 
            pointsArr[ (*numOfPoints) ] = point;
            point = createPoint( (*d) );
            (*numOfPoints)++;
            i = 0;
        }
    }
    
    if (*numOfPoints < MAX_NUMBER_OF_POINTS) {
        realloc(pointsArr, (*numOfPoints) * sizeof(Point));
    }
    freeMemPoint(point);
    fclose(input);
    return pointsArr;
}

Goal decide_command(char *arg) {
    int enumIndex;
    char *commands[] = {"spk", "wam", "ddg", "lnorm", "jacobi"};
    for (enumIndex = 0; enumIndex < MAX_CMDS; enumIndex++) {
        if (strcmp(commands[enumIndex], arg) == 0) {
            return enumIndex; /* since spk == 0 in enum Goal, the correct value will set, if found */
        } else { 
            printf("Invalid Input!\n"); /* "%s is not a goal.\nchoose from: spk / wam / ddg / lnorm / yacobi\nexits...\n", arg */ 
            assert(0);
        }
    }
    return 0;
}

Point** readPointsFromFile(int argc, char *argv[]) {
    int k, d, numOfPoints; /* max_iter? */ 
    Goal command;
    char *path;
    Point **pointsArr;
    bool findK = 0;
    
    assert( !(argc == 3) ); /* if k not provided set to 0 or exit? */ 
    k = atoi(argv[1]); 
    if (k < 0) {
        printf("Invalid Input!\n"); /* "K is not a valid integer, exits...\n" */
        assert(0);
    } else if (k == 0) {
        findK = 1;
    }
    command = decide_command(argv[2]);
    path = argv[3];
    pointsArr = readPointsArray(path, &d, &numOfPoints);
    if (k > numOfPoints) {
        printf("Invalid Input!\n"); /* "K is not smaller then n, exits...\n" */ 
        assert(0);
    }
    printPointsArr(pointsArr, numOfPoints);
    printf("command: %d\tfindK: %d\n", command, findK); /* set for gcc no to cry */ 
    return pointsArr;
}

int main(int argc, char *argv[]) {
    readPointsFromFile(argc, argv);
    return 0;
}