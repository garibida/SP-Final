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

typedef struct {
    double value;
    Point* vector;
} Eigen;

typedef struct {
    Eigen* arr;
    int length;
} Eigens_Arr;


/* ! validated, ? coded */ 
int validateInput(); /* !? */ 

/* Algorith's operations section */
Matrix* computeMatrixW(Point** pointsArr, int n);
Matrix* computeMatrixD(Matrix *W);
Matrix* computeMatrixDMinusHalf(Matrix *D);
Matrix* computeMatrixL(Matrix *W, Matrix *D); 
Matrix* computeMatrixLnorm(Matrix *L, Matrix *D); /* write tester */ 
int eigengapGetK(Eigens_Arr* eigens); /* write tester */ 

/* Point's operations section */
Point* createPoint(int d);
Point* setDataPointVal(Point *point, int index, double value);
double getDataPointVal(Point *point, int index);
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
Matrix* createUnitMatrix(int dim, bool isSymmetric);
Matrix_data createMatrixData(int rows, int cols);
Matrix_data createSymmetricMatrixData(int rows);
Matrix* cloneMatrix(Matrix* A);
void updateMatrixSymmertircStatus(Matrix* A);
double getMatrixValue(Matrix* A, int row, int col);
void setMatrixValue(Matrix* A, int row, int col, double value);
void multiply_scalar(Matrix *A, double scalar);
Matrix* add(Matrix *A, Matrix *B);
Matrix* sub(Matrix *A, Matrix *B);
Matrix* multiply(Matrix* A, Matrix* B);
bool isMatrixEqual(Matrix *A, Matrix *B);
void printMatrix(Matrix* A);
Point* createPointFromMatrixCol(Matrix* A, int col);
int compareEigens(const void *a, const void *b);
Eigens_Arr* getSortedEigen(Matrix* A);

/* Jacobi algorithm */
typedef struct
{
    int i;
    int j;
    double value;
} MaxAbsulteValue;

MaxAbsulteValue getmaxAbsulteValue(Matrix* A);
double getTheta(Matrix* A, MaxAbsulteValue mav);
double getT(double theta);
double getC(double t);
Matrix* createP(int dim, double c, double s, MaxAbsulteValue mav);
Matrix* createAtag(Matrix* A, double c, double s, MaxAbsulteValue mav);
double getOffDiagMatrixSquareSum(Matrix* A);
bool isNeedToStopJabobi(Matrix* A, Matrix* Atag);
Matrix* jacobiAlgo(Matrix** A_origin);
