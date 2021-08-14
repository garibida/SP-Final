#define true 1
#define false 0
#define TestMode 1

#define MatrixIterRows(A, i) for ((i) = 0; (i) < ((A) -> rows); (i)++)
#define MatrixIterCols(A, j) for ((j) = 0; (j) < ((A) -> cols); (j)++)
#define MatrixIterColsSym(A, i, j) for ((j) = 0; (j) <= (i); (j)++)
#define MAX_CMDS 3
#define MAX_NUMBER_OF_POINTS 1000
#define MAX_FEATURES 10
#define EPSILON 0.0001 /* set to 4 digits after the dot */ 


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
    Point** points;
    int n;
} PointsArray;

typedef struct {
    double value;
    Point* vector;
} Eigen;

typedef struct {
    Eigen* arr;
    int length;
} Eigens_Arr;

typedef enum { 
    spk = 0, 
    wam = 1, 
    ddg = 2, 
    lnorm = 3, 
    jacobi = 4
} Goal;

struct node
{
    Point* point;
    struct node *next;
};

struct linked_list
{
    struct node *head;
    struct node *tail;
    int length;
};

typedef struct node node;
typedef struct linked_list linked_list;

/* get input and validation */
Goal decide_command(char *arg);
PointsArray* readPointsArray(char *path);
PointsArray* readPointsFromFile(int argc, char *argv[]);

/* Algorith's operations section */
Matrix* computeMatrixW(PointsArray *pointsArr);
Matrix* computeMatrixD(Matrix *W);
Matrix* computeMatrixDMinusHalf(Matrix *D);
Matrix* computeMatrixL(Matrix *W, Matrix *D); 
Matrix* computeMatrixLnorm(Matrix *L, Matrix *D); 
int eigengapGetK(Eigens_Arr* eigens);
Matrix* computeMatrixU(Eigens_Arr* eigens, int k);
double* getRowsSqureRootSum(Matrix* U);
Matrix* computeMatrixT(Matrix* U);

/* Point's operations section */
Point* createPoint(int d);
void setDataPointVal(Point *point, int index, double value);
double getDataPointVal(Point *point, int index);
void printPoint(Point* point);
int isPointsEquel(Point *point1, Point* point2);
double computeDist(Point *point1, Point* point2);
double computeDistW(Point *point1, Point* point2); 
void freeMemPoint(Point *point);
Point* copy_point(Point *point);
/* Points Arr */
PointsArray* createPointsArr(int n);
Point* getPointFromArr(PointsArray* pointsArr, int i);
void setPointInArr(PointsArray* pointsArr, int i, Point* point);
void reallocPointsArr(PointsArray* pointsArr, int n);
void printPointsArr(PointsArray *pointsArr);
void freeMemPointsArr(PointsArray *pointsArr);

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
Point* createPointFromMatrixRow(Matrix* A, int row);
int compareEigens(const void *a, const void *b);
Eigens_Arr* getSortedEigen(Matrix* A);
Point* convertMatrixRowsToPoints(Matrix* A);

/* Jacobi algorithm */
typedef struct
{
    int i;
    int j;
    double value;
} MaxAbsulteValue;

MaxAbsulteValue getMaxAbsulteValue(Matrix* A);
double getTheta(Matrix* A, MaxAbsulteValue mav);
double getT(double theta);
double getC(double t);
Matrix* createP(int dim, double c, double s, MaxAbsulteValue mav);
Matrix* createAtag(Matrix* A, double c, double s, MaxAbsulteValue mav);
double getOffDiagMatrixSquareSum(Matrix* A);
bool isNeedToStopJabobi(Matrix* A, Matrix* Atag);
void calacJacobiParams(Matrix* A, MaxAbsulteValue mav, double *c, double *s);
void calcJacobiV(Matrix* A, MaxAbsulteValue mav, double c, double s, Matrix** V);
Matrix* jacobiAlgo(Matrix** A_origin);

/* K - Means */
PointsArray* kmeans(PointsArray *pointsArr, PointsArray *centroidsArr, int k, int max_iter);
bool computeCluster(int k, PointsArray *centroidsArr, PointsArray *pointsArr);
bool computeNewCentroids(linked_list** clusters, PointsArray *centroidsArr, int k);
void printOutput(double** centroids, int k, int d);
void freeDouble2DArray(double **centroids, int k);

/* Link List */
void addToList(linked_list* list, Point* point);
void freeList(linked_list* list, int isDeletePoint);
void freeNode(node* n, int isDeletePoint);