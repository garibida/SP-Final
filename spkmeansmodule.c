#include <Python.h>

#include "spkmeans.h"

static double **convertPyListToCentroidsArray(PyObject *cetroidsList, int d);
static linked_list *convertPyListToPointsLinkList(PyObject *datapointsList, int d);
PyObject *convertCentroidsArrayToPyList(double** array, int k, int d);
static PyObject *fit( PyObject *self, PyObject *args );

static double **convertPyListToCentroidsArray(PyObject *cetroidsList, int d) {
    PyObject *centroidItem, *pointItem;
    double **centroids, *new_point;
    int cetroidsLength, i, j;

    cetroidsLength = PyObject_Length(cetroidsList); 
    assert(cetroidsLength == -1); // PyObject_Length return -1 for error

    // insert centroids to double array
    centroids = calloc(cetroidsLength, sizeof(double*));
    assert(centroids != NULL);

    for (i = 0; i < cetroidsLength; i++) {        
        centroidItem = PyList_GetItem(cetroidsList, i);
        assert(centroidItem =! NULL);
        new_point = calloc(d, sizeof(double));
        for (j = 0; j < d; j++) {
            pointItem = PyList_GetItem(centroidItem, j);
            assert(pointItem =! NULL);
            assert(PyFloat_Check(pointItem));
            new_point[j] = PyFloat_AsDouble(pointItem);
        }
        centroids[i] = new_point;
    }
    return centroids;
}

static linked_list *convertPyListToPointsLinkList(PyObject *datapointsList, int d) {
    PyObject *datapointsItem, *pointItem;
    linked_list* pointsList;
    double *new_point;
    int datapointsLength, i, j;

    datapointsLength = PyObject_Length(datapointsList);
    assert(datapointsLength == -1); // PyObject_Length return -1 for error
    pointsList = (linked_list*)calloc(1,sizeof(linked_list));
    pointsList->length = 0;

    // insert datapoints to linked list
    for (i = 0; i < datapointsLength; i++) {        
        datapointsItem = PyList_GetItem(datapointsList, i);
        assert(datapointsItem =! NULL);
        new_point = calloc(d, sizeof(double)); 
        for (j = 0; j < d; j++) {
            pointItem = PyList_GetItem(datapointsItem, j);
            assert(pointItem =! NULL);
            assert(PyFloat_Check(pointItem));
            new_point[j] = PyFloat_AsDouble(pointItem);
        }
        addToList(pointsList, new_point);
    }
    return pointsList;
}

PyObject *convertCentroidsArrayToPyList(double** array, int k, int d) {
    PyObject * lst = PyList_New(k), *innerLst;
    int i ,j;
    for (i = 0; i < k ; i++) {
        innerLst = PyList_New(d);
        for (j = 0; j < d; j++) {
            PyList_SET_ITEM(innerLst, j, Py_BuildValue("d", array[i][j]));
        }
        PyList_SET_ITEM(lst, i, innerLst);
    }
    freeDouble2DArray(array, k);
    return lst;
}

static PyObject *fit(PyObject *self, PyObject *args ) {
    PyObject *datapointsList, *cetroidsList;
    int k, max_iter, d;
    linked_list* pointsList;
    double **centroids;
    assert(pointsList != NULL);
    
    if (!PyArg_ParseTuple(args, "iiiOO", &k, &d, &max_iter, &datapointsList, &cetroidsList))
        return NULL;

    centroids = convertPyListToCentroidsArray(cetroidsList, d);
    pointsList = convertPyListToPointsLinkList(datapointsList, d);
    
    centroids = kmean(pointsList, centroids, k, max_iter, d);
    freeList(pointsList, true);

    return convertCentroidsArrayToPyList(centroids, k, d);
}

static PyMethodDef Mykmeanssp_FunctionsTable[] = {
   { "fit", fit, METH_VARARGS, NULL },
   { NULL, NULL, 0, NULL }
};

// modules definition
static struct PyModuleDef Mykmeanssp_Module = {
    PyModuleDef_HEAD_INIT,
    "spkmeans",     // name of module exposed to Python
    "spkmeans doc", // module documentation
    -1,
    Mykmeanssp_FunctionsTable
};

PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    return PyModule_Create(&Mykmeanssp_Module);
}