#include "spkmeans.c"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h> /* check */

/* test section */
void testMain(bool isDebug);
void testMultiplyMatrixs(bool isDebug);
Point** pointsForTestE0();
void TestE0(bool isDebug);
Point** pointsForTestE1();
void TestE1(bool isDebug);
void testJacobi(bool isDebug);
void testEigen(bool isDebug);

# define isDoubleEqual(_a, _b) (fabs((_a) - (_b)) < 0.0001)


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
    Matrix *A , *U;
    bool testResult;
    double epsilon = 0.0000001;
    Eigens_Arr *eigensArr;
    Eigen *eigens, *expectedEigens1, *expectedEigens2;
    int i, j, k;
    A = createMatrix(3, 3, true);
    setMatrixValue(A, 0, 0, 3.0);
    setMatrixValue(A, 1, 0, 2.0);
    setMatrixValue(A, 1, 1, 0.0);
    setMatrixValue(A, 2, 0, 4.0);
    setMatrixValue(A, 2, 1, 2.0);
    setMatrixValue(A, 2, 2, 3.0);

    expectedEigens1 = (Eigen *) calloc(3, sizeof(Eigen));
    assert(expectedEigens1 != NULL);
    expectedEigens1[0].value = -1;
    expectedEigens1[0].vector = createPoint(3);
    expectedEigens1[1].value = -1;
    expectedEigens1[1].vector = createPoint(3);
    expectedEigens1[2].value = 8;
    expectedEigens1[2].vector = createPoint(3);
    setDataPointVal(expectedEigens1[0].vector, 0, 1 / sqrt(2));
    setDataPointVal(expectedEigens1[0].vector, 1, 0);
    setDataPointVal(expectedEigens1[0].vector, 2, - 1 / sqrt(2));
    setDataPointVal(expectedEigens1[1].vector, 0, - 1 / (3 * sqrt(2)));
    setDataPointVal(expectedEigens1[1].vector, 1, 4 / (3 * sqrt(2)));
    setDataPointVal(expectedEigens1[1].vector, 2, - 1 / (3 * sqrt(2)));
    setDataPointVal(expectedEigens1[2].vector, 0, 2.0 / 3.0);
    setDataPointVal(expectedEigens1[2].vector, 1, 1.0 / 3.0);
    setDataPointVal(expectedEigens1[2].vector, 2, 2.0 / 3.0);

    expectedEigens2 = (Eigen *) calloc(3, sizeof(Eigen));
    assert(expectedEigens2 != NULL);
    expectedEigens2[0].value = -1;
    expectedEigens2[0].vector = createPoint(3);
    expectedEigens2[1].value = -1;
    expectedEigens2[1].vector = createPoint(3);
    expectedEigens2[2].value = 8;
    expectedEigens2[2].vector = createPoint(3);
    setDataPointVal(expectedEigens2[1].vector, 0, 1 / sqrt(2));
    setDataPointVal(expectedEigens2[1].vector, 1, 0);
    setDataPointVal(expectedEigens2[1].vector, 2, - 1 / sqrt(2));
    setDataPointVal(expectedEigens2[0].vector, 0, - 1 / (3 * sqrt(2)));
    setDataPointVal(expectedEigens2[0].vector, 1, 4 / (3 * sqrt(2)));
    setDataPointVal(expectedEigens2[0].vector, 2, - 1 / (3 * sqrt(2)));
    setDataPointVal(expectedEigens2[2].vector, 0, 2.0 / 3.0);
    setDataPointVal(expectedEigens2[2].vector, 1, 1.0 / 3.0);
    setDataPointVal(expectedEigens2[2].vector, 2, 2.0 / 3.0);

    eigensArr = getSortedEigen(A);
    eigens = eigensArr->arr;
    testResult = true;
    for(i = 0; i < 3; i++) {
        if (fabs(eigens[i].value - expectedEigens1[i].value) > epsilon) {
            testResult = false;
            break;
        }
        if(!isPointsEquel(eigens[i].vector, expectedEigens1[i].vector) &&
           !isPointsEquel(eigens[i].vector, expectedEigens2[i].vector)) {
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

    (eigengapGetK(eigensArr) == 1) ?
        printf("'test Eigengap Heuristic'\tresult: Great!\n") :
        printf("'test Eigengap Heuristic'\tresult: Problem!\n");

    k = 1;
    U = computeMatrixU(eigensArr, k);

    assert(U->cols = k);
    assert(U->rows = 3);
    testResult = true;
    MatrixIterRows(U, i) {
        MatrixIterCols(U, j) {
            if (!isDoubleEqual(getMatrixValue(U, i, j), getDataPointVal((eigensArr->arr)[j].vector, i))) {
                printf("im\n");
                testResult = false;
                break;
            }
        }
    }
    (testResult) ?
        printf("'test Matrix U'\t\t\tresult: Great!\n") :
        printf("'test Matrix U'\t\t\tresult: Problem!\n");
}

void testMain(bool isDebug) {
    testMultiplyMatrixs(isDebug);
    testJacobi(isDebug);
    testEigen(isDebug);
    TestE0(isDebug);
    TestE1(isDebug);
}

int main() {
    if (TestMode) {
        testMain(false);
    }
    return 0;
}