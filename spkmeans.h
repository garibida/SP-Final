
/* ######################################################################### */
/*                                   MACROs                                  */
/* ######################################################################### */
#define MatrixIterRows(A, i)                    for ((i) = 0; (i) < ((A) -> rows); (i)++)
#define MatrixIterCols(A, j)                    for ((j) = 0; (j) < ((A) -> cols); (j)++)
#define MatrixIterColsSym(A, i, j)              for ((j) = 0; (j) <= (i); (j)++)
#define MatrixIterColsSymUperTriengle(A, i, j)  for ((j) = (i); (j) < ((A) -> cols); (j)++)
#define ASSERT_M(cond, msg)                     if ( !(cond) ) {printf(msg); assert(false);}

/* ######################################################################### */
/*                                   Consts                                  */
/* ######################################################################### */

#define MAX_CMDS                                3
#define ENUM_COUNT                              5
#define MAX_NUMBER_OF_POINTS                    1000 /* set to 50 ######################################################################################### TODO: change */ 
#define MAX_FEATURES                            1000 /* set to 10 ######################################################################################### TODO: change */ 
#define EPSILON                                 0.0001 /* set to 4 digits after dot */ 
#define EPSILON_YACOBI                          1.0e-15
#define true                                    1
#define false                                   0
#define ERROR_MSG                               "An Error Has Occured\n"
#define INVALID_INPUT_MSG                       "Invalid Input!\n"

/* ######################################################################### */
/*                             Structs and Enums                             */
/* ######################################################################### */

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
    int index;
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

typedef struct
{
    int i;
    int j;
    double value;
} MaxAbsulteValue;

/* ######################################################################### */
/*                                 Functions                                 */
/* ######################################################################### */

/* Matrix's operations section */
Matrix* createMatrix(int rows, int cols, bool isSymmetric);
Matrix_data createMatrixData(int rows, int cols);
Matrix_data createSymmetricMatrixData(int rows);
Matrix* createUnitMatrix(int dim, bool isSymmetric);
Matrix* cloneMatrix(Matrix* A);
void updateMatrixSymmertircStatus(Matrix* A);
double getMatrixValue(Matrix* A, int row, int col);
void setMatrixValue(Matrix* A, int row, int col, double value);
Matrix* sub(Matrix *A, Matrix *B, bool doFree);
Matrix* multiply(Matrix* A, Matrix* B, bool doFree);
void freeMatrix(Matrix* A);
void printMatrix(Matrix* A);
Point* createPointFromMatrixCol(Matrix* A, int col);
Point* createPointFromMatrixRow(Matrix* A, int row);
PointsArray* matrixToPointsArray(Matrix *A);
Matrix* PointsArrayToMatrix(PointsArray *pointsArr);

/* Eigen's operations section */
int compareEigens(const void *a, const void *b);
Eigens_Arr* getEigens(Matrix **A);
Eigens_Arr* getSortedEigens(Matrix **A);
void freeEigens(Eigens_Arr *eigens);
void printEigens(Eigens_Arr *eigens);

/* Point's operations section */
Point* createPoint(int d);
Point* createPointWithVals(int d, double *values);
void setDataPointVal(Point *point, int index, double value);
double getDataPointVal(Point *point, int index);
void printPoint(Point *point, bool isLast);
int isPointsEquel(Point *point1, Point* point2);
void freeMemPoint(Point *point);
Point* copy_point(Point *point);
double computeDist(Point *point1, Point* point2);

/* Points Array's operations section */
PointsArray* createPointsArr(int n);
Point* getPointFromArr(PointsArray* pointsArr, int i);
void setPointInArr(PointsArray* pointsArr, int i, Point* point);
void reallocPointsArr(PointsArray* pointsArr, int n);
void printPointsArr(PointsArray *pointsArr);
void freeMemPointsArr(PointsArray *pointsArr);

/* Link List */
void addToList(linked_list* list, Point* point);
void freeList(linked_list* list, int isDeletePoint);
void freeNode(node* n, int isDeletePoint);

/* Algorith's operations section */
double computeDistW(Point *point1, Point* point2); 
Matrix* computeMatrixW(PointsArray *pointsArr);
Matrix* computeMatrixD(Matrix *W);
Matrix* computeMatrixDMinusHalf(Matrix *D);
Matrix* computeMatrixLnorm(Matrix *L, Matrix *D); 
int eigengapGetK(Eigens_Arr* eigens);
Matrix* computeMatrixU(Eigens_Arr* eigens, int k);
double* getRowsSqureRootSum(Matrix* U);
Matrix* computeMatrixT(Matrix* U);

/* Jacobi algorithm */
MaxAbsulteValue getMaxAbsulteValue(Matrix* A);
double getTheta(Matrix* A, MaxAbsulteValue mav);
double getT(double theta);
double getC(double t);
Matrix* createP(int dim, double c, double s, MaxAbsulteValue mav);
Matrix* createAtag(Matrix* A, double c, double s, MaxAbsulteValue mav);
double getOffDiagMatrixSquareSum(Matrix* A);
bool isNeedToStopJacobi(Matrix* A, Matrix* Atag);
void calacJacobiParams(Matrix* A, MaxAbsulteValue mav, double *c, double *s);
void calcJacobiV(Matrix* A, MaxAbsulteValue mav, double c, double s, Matrix** V);
Matrix* jacobiAlgo(Matrix** A_origin);

/* K - Means */
PointsArray* kmeans(PointsArray *pointsArr, PointsArray *centroidsArr, int k, int max_iter);
bool computeCluster(int k, PointsArray *centroidsArr, PointsArray *pointsArr);
bool computeNewCentroids(linked_list** clusters, PointsArray *centroidsArr, int k);
PointsArray* getIntialCentroids(PointsArray *pointsArr, int k);

/* Main Functions */
void matrixPrinter(PointsArray *points ,Goal goal);
int doSpk(PointsArray **points, int k);
PointsArray* readPointsArray(char *path);
Goal decide_command(char *arg);