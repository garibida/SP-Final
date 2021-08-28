#include <Python.h>
#include "spkmeans.h"

static PointsArray *convertPyListToPointsArray(PyObject *datapointsList);
PyObject *convertPointsArrayToPyList(PointsArray *pointsArray);
PyObject *convertCentroidsArrayToPyList(double** array, int k, int d);
static PyObject *fit(PyObject *self, PyObject *args);
static PyObject *doSpkPython(PyObject *self, PyObject *args);
static PyObject *printMatrixes(PyObject *self, PyObject *args);

static PointsArray *convertPyListToPointsArray(PyObject *datapointsList) {
    PyObject *datapointsItem, *pointItem;

    PointsArray* pointsArr;
    Point *newPoint;

    int n, d, i, j;

    n = PyObject_Length(datapointsList);
    assert(n == -1); // PyObject_Length return -1 for error

    pointsArr = createPointsArr(n);

    // insert datapoints to linked list
    for (i = 0; i < n; i++) {
        datapointsItem = PyList_GetItem(datapointsList, i);
        assert(datapointsItem =! NULL);

        d = PyObject_Length(datapointsItem);
        assert(d == -1); // PyObject_Length return -1 for error

        newPoint = createPoint(d);
        for (j = 0; j < d; j++) {
            pointItem = PyList_GetItem(datapointsItem, j);
            assert(pointItem =! NULL);
            assert(PyFloat_Check(pointItem));
            setDataPointVal(newPoint, j, PyFloat_AsDouble(pointItem));
        }
        setPointInArr(pointsArr, i, newPoint);
    }
    return pointsArr;
}

PyObject *convertPointsArrayToPyList(PointsArray *pointsArray) {
    PyObject *lst, *innerLst;
    Point *point;
    lst = PyList_New(pointsArray->n);
    int i ,j;
    for (i = 0; i < (pointsArray->n); i++) {
        point = getPointFromArr(pointsArray, i);
        innerLst = PyList_New(point->d);
        for (j = 0; j < (point->d); j++) {
            PyList_SET_ITEM(innerLst, j, Py_BuildValue("d", getDataPointVal(point, j)));
        }
        PyList_SET_ITEM(lst, i, innerLst);
    }
    freeMemPointsArr(pointsArray);
    return lst;
}

static PyObject *fit(PyObject *self, PyObject *args) {
    PyObject *datapointsList, *cetroidsList;
    int k, max_iter = 300;
    PointsArray *points, *centroids;

    if (!PyArg_ParseTuple(args, "iOO", &k, &datapointsList, &cetroidsList))
        return NULL;

    centroids = convertPyListToPointsArray(cetroidsList);
    points = convertPyListToPointsArray(datapointsList);

    centroids = kmeans(points, centroids, k, max_iter);
    freeMemPointsArr(points);
    printCentroids(centroids);
    freeMemPointsArr(centroids);

    Py_RETURN_NONE;
}

static PyObject *doSpkPython(PyObject *self, PyObject *args) {
    PyObject *datapointsList, *lst;
    int k;
    PointsArray *points;
    assert(pointsList != NULL); /* what is pointsList? */
    
    if (!PyArg_ParseTuple(args, "iO", &k, &datapointsList))
        return NULL;

    points = convertPyListToPointsArray(datapointsList);

    k = doSpk(&points, k);
    datapointsList = convertPointsArrayToPyList(points);

    lst = PyList_New(2);
    PyList_SET_ITEM(lst, 0, Py_BuildValue("i", k));
    PyList_SET_ITEM(lst, 1, Py_BuildValue("O", datapointsList));

    return lst;
}

static PyObject *printMatrixes(PyObject *self, PyObject *args) {
    PyObject *datapointsList;
    Goal goal;
    int goal_int;
    PointsArray *points;
    
    if (!PyArg_ParseTuple(args, "iO", &goal_int, &datapointsList))
        return NULL;

    goal = (Goal) goal_int;
    points = convertPyListToPointsArray(datapointsList);

    matrixPrinter(points, goal);

    Py_RETURN_NONE;
}

static PyMethodDef SpkmeansModule_FunctionsTable[] = {
   { "fit", fit, METH_VARARGS, NULL },
   { "spk", doSpkPython, METH_VARARGS, NULL },
   { "printMatrixes", printMatrixes, METH_VARARGS, NULL },
   { NULL, NULL, 0, NULL }
};

// modules definition
static struct PyModuleDef SpkmeansModule = {
    PyModuleDef_HEAD_INIT,
    "spkmeansModule",     // name of module exposed to Python
    "spkmeansModule doc", // module documentation
    -1,
    SpkmeansModule_FunctionsTable
};

PyMODINIT_FUNC PyInit_spkmeansModule(void) {
    return PyModule_Create(&SpkmeansModule);
}