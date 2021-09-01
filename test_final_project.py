import math
import os
import pandas as pd
import numpy as np
from math import e
import os.path
import filecmp
import time

epsilon = 1e-6
C_file_name = "hw1"


def get_vectors(filename):
    result = np.genfromtxt(filename, delimiter=',').tolist()
    if type(result) == float:
        result = [[result]]
    elif type(result[0]) == float:
        result = [result]
    return result


def build_weight_matrix(vectors):
    N_size = vectors.shape[0]
    result = np.fromfunction(lambda i, j: e ** (calc_norm(vectors[i], vectors[j])/(-2)), (N_size, N_size), dtype=int)
    return result - np.eye(N_size)


def build_diagonal_matrix(vectors):
    wMatrix = build_weight_matrix(vectors)
    N_size = wMatrix.shape[0]
    result = np.zeros((N_size, N_size))
    for i in range(len(wMatrix)):
        result[i][i] = np.sum(wMatrix[i])
    return result


def build_lnorm_matrix(vectors):
    wMatrix = build_weight_matrix(vectors)
    dMatrix = build_diagonal_matrix(vectors)
    for i in range(len(dMatrix)):
        dMatrix[i][i] = 1 / np.sqrt(dMatrix[i][i])
    return np.eye(dMatrix.shape[0]) - dMatrix @ wMatrix @ dMatrix


# def build_jacobi_matrix(lmatrix):
#     epsilon = pow(10,-15)
#     jacobi_values = lmatrix
#     n_size = len(jacobi_values)
#     jacobi_vectors = np.eye(n_size)
#     cur_f_norm = calc_off(jacobi_values)
#     for iter in range(100):
#         i, j = find_off_diag_max(jacobi_values)
#         psi = (jacobi_values[j][j] - jacobi_values[i][i]) / (2 * jacobi_values[i][j])
#         t = (-1 if psi < 0 else 1) / (abs(psi) + math.sqrt(psi ** 2 + 1))
#         c = 1 / math.sqrt(1 + t ** 2)
#         s = t * c
#         P_matrix = build_P_matrix(n_size,i,j,c,s)
#         jacobi_values = P_matrix.transpose() @ jacobi_values @ P_matrix
#         jacobi_vectors = jacobi_vectors @ P_matrix
#         prev_f_norm = cur_f_norm
#         cur_f_norm = calc_off(jacobi_values)
#         if abs(cur_f_norm - prev_f_norm) < epsilon:
#             break
#     return jacobi_values, jacobi_vectors


# def build_P_matrix(N, i, j, c, s):
#     result = np.eye(N)
#     result[i][i] = c
#     result[j][j] = c
#     result[i][j] = s
#     result[j][i] = -s
#     return result


# def calc_off(matrix):
#     result = 0
#     for i in range(len(matrix)):
#         for j in range(i+1, len(matrix)):
#             result += matrix[i][j]**2 + matrix[j][i]**2
#     return result


# def find_off_diag_max(matrix):
#     matrix_max = -1
#     max_i = -1
#     max_j = -1
#     for i in range(len(matrix)):
#         for j in range(i+1, len(matrix)):
#             if abs(matrix[i][j]) > matrix_max:
#                 matrix_max = abs(matrix[i][j])
#                 max_i = i
#                 max_j = j
#     return max_i, max_j


def calc_norm(x1, x2):
    return np.sqrt(np.sum((x1 - x2)**2, axis=2))



def print_list(lst):
    for i in range(len(lst) - 1):
        print("{:.4f}".format(lst[i]), end=",")
    print("{:.4f}".format(lst[len(lst) - 1]))

def print_mat(mat):
    for line in mat:
        print_list(line)




def generate_correct_files_general():
    print("======Generating correct files genreal=========")
    print("======SPK WILL BE ALWAYS CORRECT=========")
    exec = ["python3 spkmeans.py", os.path.join(".", C_file_name)]
    filesToCreate = pd.read_csv(os.path.join(".", "tests", "FilesToCreate.csv"))
    input_index = 0
    goals = ["spk", "ddg", "wam", "lnorm"]
    for row in filesToCreate.itertuples():
        k_orig = row.n_clusters
        input_data_filename = os.path.join(".", "tests", "test_data", "spk_tests", f"test{input_index}.csv")     
        for goal in goals: 
            if goal == "spk":
                k_values = [k_orig, k_orig//2, 0]
            else:
                k_values = [1]       
            for k in k_values:
                for ex in exec:
                    lng = "P" if ex == "python3 spkmeans.py" else "C"
                    result_file = f"test{input_index}_{goal}_{k}_output_{lng}.txt"
                    result_path = os.path.join(".", "tests", "reference", "general", result_file)
                    start = time.time()
                    returnCode = os.system(f"{ex} {k} {goal} {input_data_filename} > {result_path}")
                    end = time.time()
                    if returnCode != 0:
                        print("Something went wrong generating this file...")
                        print(f"k = {k} goal = {goal} lng = {lng} input = {input_index}")
        input_index+=1

def generate_correct_files_jacobi():
    print("======Generating correct files jacobi=========")
    print("======JACOBI WILL BE ALWWAYS TRUE=========")
    exec = ["python3 spkmeans.py", os.path.join(".", C_file_name)]
    
    for i in range(11):
        input_data_filename = os.path.join(".", "tests", "test_data", "jacobi_tests", f"test{i}.csv")     
        for ex in exec:
            lng = "P" if ex == "python3 spkmeans.py" else "C"
            result_file = f"test{i}_jacobi_output_{lng}.txt"
            result_path = os.path.join(".", "tests", "reference", "jacobi", result_file)
            start = time.time()
            returnCode = os.system(f"{ex} 1 jacobi {input_data_filename} > {result_path}")
            end = time.time()
            if returnCode != 0:
                print("Something went wrong generating this file...")
                print(f"k = 1 goal = jacobi lng = {lng} input = {input_index}")

                
               
def handle_goal(ex, k_values, goal, input_data_filename, result_file, input_index):
    errors = 0
    fails_numeric = 0 
    failes_file = 0
    success = 0
    lng = "P" if ex == "python3 spkmeans.py" else "C"
    for k in k_values:
        start = time.time()
        returnStatus = os.system(f"{ex} {k} {goal} {input_data_filename} > {result_file}")
        end = time.time()
        result_matrix = get_vectors(result_file)
        vectors = np.array(get_vectors(input_data_filename))
        correctOutputPath = os.path.join(".", "tests", "reference", "general", f"test{input_index}_{goal}_{k}_output_{lng}.txt")
        if returnStatus != 0:
            errors+=1
        if goal == "spk":
            correct_matrix = get_vectors(correctOutputPath)
        elif goal == "wam":
            correct_matrix = build_weight_matrix(vectors)
        elif goal == "ddg":
            correct_matrix = build_diagonal_matrix(vectors)
        elif goal == "lnorm":
            correct_matrix = build_lnorm_matrix(vectors)
        correct_matrix = np.round(correct_matrix, 4).tolist()
        res_numeric = np.array_equal(correct_matrix, result_matrix)
        res_files = filecmp.cmp(result_file, correctOutputPath, shallow=False)
        if res_numeric:
            if res_files:
                res_str = "Paseed"
                success+=1
            else:
                res_str = "Passed numeric Failed file"
                failes_file+=1
        else:
            fails_numeric+=1
            if res_files:
                res_str = "Passed file fail numeric (?)"
            else:
                res_str = "Failed numeric and file"
                failes_file+=1
        time_str = "{:.4f}".format(end-start)
        print(f"check file={input_data_filename} \tlanguage={lng} \tk_value={k} \tgoal={goal} \ttime={time_str} \tResult={res_str}")
        if not res_numeric:
            print("Your output:")
            print_mat(result_matrix)
            print("correct output:")
            print_mat(correct_matrix)
            
    return errors, fails_numeric, failes_file, success

if __name__ == "__main__":
    print("======Generating input from FilesToCreate.csv and general jacobi files=============") 
    exec(open(os.path.join(".", "tests", "CreateTestFiles.py")).read())
    testFiles = pd.read_csv(os.path.join(".", "tests", "FilesToCreate.csv"))
    print(f"created {testFiles.shape[0]} general test files")
    print("created 11 jacobi test files")
    
    # uncomment this next two line to generate your files as reference for correct files
    # generate_correct_files_general()
    # generate_correct_files_jacobi()

    print("===================Comparing results======================")
    exec = ["python3 spkmeans.py", os.path.join(".", C_file_name)]
    result_file = "tmp.txt"
    correctFiles = 0
    wrong_numeric = 0
    wrong_files = 0
    erroeFiles = 0

    filesToCreate = pd.read_csv(os.path.join(".", "tests", "FilesToCreate.csv"))
    input_index = 0

    for row in filesToCreate.itertuples():
        k_orig = row.n_clusters
        goals = ["spk", "ddg", "wam", "lnorm"]
        input_data_filename = os.path.join(".", "tests", "test_data", "spk_tests", f"test{input_index}.csv")
        for ex in exec:
            for goal in goals:
                if goal == "spk":
                    k_values = [k_orig, k_orig//2, 0]
                else:
                    k_values = [1]       
                result = handle_goal(ex, k_values, goal, input_data_filename, result_file, input_index)
                erroeFiles+=result[0]
                wrong_numeric+=result[1]
                wrong_files+=result[2]
                correctFiles+=result[3]
        input_index+=1

    # Jacobi testing
    for i in range(11):
        input_data_filename = os.path.join(".", "tests", "test_data", "jacobi_tests", f"test{i}.csv")
        for ex in exec:
            lng = "P" if ex == "python3 spkmeans.py" else "C"
            start = time.time()
            returnStatus = os.system(f"{ex} 1 jacobi {input_data_filename} > {result_file}")
            end = time.time()
            result_matrix = get_vectors(result_file)
            if returnStatus != 0:
                erroeFiles+=1
                break
            correct_matrix = get_vectors(os.path.join(".", "tests", "reference", "jacobi", f"test{i}_jacobi_output_{lng}.txt"))
            res = np.array_equal(correct_matrix, result_matrix)
            res_str = "Passed" if res else "Failed"
            time_str = "{:.4f}".format(end-start)
            print(f"check file={input_data_filename} \tlanguage={lng} \tgoal=jacobi \ttime={time_str} \tResult={res_str}")
            if not res:
                print("Your output:")
                print_mat(result_matrix)
                print("correct output:")
                print_mat(correct_matrix)
                wrong_numeric+=1
            else:
                correctFiles+=1

 
    print("=======================SUMMARY==============================")
    print(f"Total result:\n")
    print(f"\tTests passed: {correctFiles}\n")
    print(f"\tTests failed numeric: {wrong_numeric}\n")
    print(f"\tTests failed text: {wrong_files}\n")
    print(f"\tTests errored: {erroeFiles}\n")

        

