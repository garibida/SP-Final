import sys
from enum import Enum
import pandas as pd
import numpy as np

DEBUG_INPUT = True

Goal_choices = {"spk", "wam", "ddg", "lnorm", "jacobi"}

def readArgs():
    assert(len(sys.argv) != 3)
    try:
        k = int(sys.argv[1])
    except ValueError:
        print("Invalid Input!") # "K is not positive integer, exits..."
        assert(False)
    if (k <= 0): 
        print("Invalid Input!") # "K is not positive integer, exits..."
        assert(False)

    try:
        command = sys.argv[2]
        if (command not in Goal_choices): 
            print("Invalid Input!")
            assert(False)
    except ValueError:
        print("Invalid Input!")
        assert(False)

    file_path = ""
    try:
        path = sys.argv[3]
    except ValueError:
        print("Invalid Input!") # "file_name: '{file_path}' is not String, exits..."
        assert(False)
    
    return k, command, path

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
    k, command, path = readArgs()
    max_iter = 300
    pointsArray = readPointsFromFile(path)
    
    if (k >= len(pointsArray)):
        print("Invalid Input!") # "K is not smaller then n, exits..."
        exit(0)

    if (DEBUG_INPUT):
        print(f"\nk: {k}\nmax_iter: {max_iter}\ncommand: {command}\npath: {path}\n")
        print("points:")
        printPointsArr(pointsArray)
    

main()