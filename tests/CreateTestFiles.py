import pandas as pd
import numpy as np
from sklearn.datasets import make_blobs
import matplotlib.pyplot as plt
import os


'''Creates the files for all goals except jacobi'''

filesToCreate = pd.read_csv(os.path.join(".", "tests", "FilesToCreate.csv"))
i = 0

for row in filesToCreate.itertuples():
    if (i == 10):
        continue
    n_centers = row.n_clusters
    samples = row.n_samples
    features = row.features
    
    file_name = f"test{i}.csv"
    
    X, y = make_blobs(n_samples=samples, centers=n_centers, n_features=features,
              random_state=31, shuffle=True, cluster_std=0.3)
    resultPath = os.path.join(".", "tests", "test_data", "spk_tests", file_name)
    np.savetxt(resultPath, X, delimiter=",", fmt="%.4f")
    
    #plt.plot(X[:, 0], X[:, 1], 'o', label = 'data')
    #plt.show()
    
    i += 1
    
'''Creates the files for jacobi'''
np.random.seed(0)
dim = 2
for i in range(12):
    if (i == 11):
        continue
    mat = np.random.rand(dim, dim)
    mat = np.tril(mat) + np.tril(mat, -1).T
    resultPath = os.path.join(".", "tests", "test_data", 'jacobi_tests',f"test{i}.csv")
    np.savetxt(resultPath, mat, delimiter=",", fmt="%.4f")
    
    dim += 1
    
    