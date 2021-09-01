import sys
from enum import Enum
import pandas as pd
import numpy as np
import spkmeansModule

DEBUG_INPUT = False
MAX_ITER = 300
GOAL_CHOICES = {"spk", "wam", "ddg", "lnorm", "jacobi"}
ERROR_MSG = "An Error Has Occured\n"
INVALID_INPUT_MSG = "Invalid Input!\n"


class Goal(Enum):
    spk = 0
    wam = 1
    ddg = 2
    lnorm = 3
    jacobi = 4

def readArgs():
    assert(len(sys.argv) != 3)
    try:
        k = int(sys.argv[1])
    except ValueError:
        print(INVALID_INPUT_MSG)
        assert(False)
    if (k < 0): 
        print(INVALID_INPUT_MSG)
        assert(False)

    try:
        goal = Goal[sys.argv[2]]
    except ValueError:
        print(INVALID_INPUT_MSG)
        assert(False)

    try:
        path = sys.argv[3]
    except ValueError:
        print(INVALID_INPUT_MSG)
        assert(False)
    
    return k, goal, path

def readPointsFromFile(file_path): 
    pointsDf = pd.read_csv(file_path, header=None)
    pointsNd = pointsDf.to_numpy()
    return list(map(lambda x: x.tolist(), pointsNd))

def printPointsArr(pointsArr): 
    str = ""
    for i in range(len(pointsArr)): 
        for d in range(len(pointsArr[0])): 
            str += "{}".format(format(pointsArr[i][d], '.4f'))
            str += ","
        str = str[:-1] + "\n"
    
    print(str[:-1])

def kmeans_pp(datapoints, k):
    np.random.seed(0)
    r = np.random.choice(len(datapoints))
    centroids = [(r, datapoints[r])]
    D = [np.inf for i in range(len(datapoints))]

    Z = 1
    while(Z < k):
        for i in range(len(datapoints)):
            x = datapoints[i]
            curDist = np.linalg.norm(x - centroids[-1][1]) ** 2
            D[i] = curDist if (curDist < D[i]) else D[i]

        Z += 1
        dSum = sum(D)
        NormalizedD = list(map(lambda d: d / dSum, D))
        r = np.random.choice(len(datapoints), p=NormalizedD)
        centroids.append((r, datapoints[r]))

    return centroids

def main():
    k, goal, path = readArgs()
    max_iter = MAX_ITER # reminder if needed
    pointsArray = readPointsFromFile(path)
    
    if (k >= len(pointsArray)):
        print(INVALID_INPUT_MSG)
        assert(False)

    if (DEBUG_INPUT):
        print(f"\nk: {k}\nmax_iter: {max_iter}\ncommand: {goal}\npath: {path}\n")
        print("points:")
        printPointsArr(pointsArray)
    
    if goal is not Goal.spk:
        spkmeansModule.printMatrixes(goal.value, pointsArray)
        return

    res = spkmeansModule.spk(k, pointsArray)
    k = res[0]
    pointsArray = res[1]

    centroids = kmeans_pp(np.array(pointsArray), k)
    centroidsArray = list(map(lambda x: x[1].tolist(), centroids))
    print(",".join([str(int(c[0])) for c in centroids]))
    spkmeansModule.fit(k, pointsArray, centroidsArray)

main()