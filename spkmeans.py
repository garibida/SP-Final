import sys
from enum import Enum
import pandas as pd
import numpy as np

DEBUG_INPUT = True
MAX_ITER = 300
GOAL_CHOICES = {"spk", "wam", "ddg", "lnorm", "jacobi"}

def readArgs():
    assert(len(sys.argv) != 3)
    try:
        k = int(sys.argv[1])
    except ValueError:
        print("Invalid Input!\n") # "K is not positive integer, exits..."
        assert(False)
    if (k < 0): 
        print("Invalid Input!\n") # "K is not positive integer, exits..."
        assert(False)
    
    findK = 1 if (k == 0) else 0

    try:
        command = sys.argv[2]
        if (command not in GOAL_CHOICES): 
            print("Invalid Input!\n")
            assert(False)
    except ValueError:
        print("Invalid Input!\n")
        assert(False)

    file_path = ""
    try:
        path = sys.argv[3]
    except ValueError:
        print("Invalid Input!\n") # "file_name: '{file_path}' is not String, exits..."
        assert(False)
    
    return k, command, path, findK

def readPointsFromFile(file_path): 
    pointsDf = pd.read_csv(file_path, header=None)
    return pointsDf.sort_index().to_numpy() # sort needed?

def printPointsArr(pointsArr): 
    str = ""
    for i in range(len(pointsArr)): 
        for d in range(len(pointsArr[0])): 
            str += "{}".format(format(pointsArr[i][d], '.4f'))
            str += ","
        str = str[:-1] + "\n"
    
    print(str[:-1])

def main():
    k, command, path, findK = readArgs()
    max_iter = MAX_ITER # reminder if needed
    pointsArray = readPointsFromFile(path)
    
    if (k > len(pointsArray)):
        print("Invalid Input!") # "K is not smaller then n, exits..."
        assert(False)

    if (DEBUG_INPUT):
        print(f"\nk: {k}\nmax_iter: {max_iter}\ncommand: {command}\npath: {path}\n")
        print("points:")
        printPointsArr(pointsArray)
    
main()