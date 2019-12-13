#!/usr/bin/python3


import smacof_numerical
import numpy as np
import math

def find_hessian(weights):
    v = np.zeros(weights.shape)
    l = weights.shape[0]
    # loop over indices n < m of the (square) array weights
    for n, m in np.ndindex(weights.shape):
        if n < m:
            v[n, n] += weights[n, m]
            v[n, m] -= weights[n, m]
            v[m, m] += weights[n, m]
            v[m, n] -= weights[n, m]
    return v


def find_vplus(weights):
    hessian = find_hessian(weights)
    print("Calculating Moore-Penrose Pseudo-Inverse of the Hessian...", end = '', flush = True)
    vplus = np.linalg.pinv(hessian)
    print("...done.", flush = True)
    return vplus


def run_smacof(weights, lengths, target_dimension = 3, max_iterations = 1000):
    # subroutine smacof_embed(numvertices, target_dimension, weights, lengths, Vplus, maximum_iterations, X)
    Vplus = find_vplus(weights)
    print("Running the Fortran SMACOF algorithm...", end = '', flush = True)
    X = smacof_numerical.smacof_embed(target_dimension = target_dimension, weights = weights, lengths = lengths, Vplus = Vplus, maximum_iterations = maximum_iterations)
    print("...done.")

    return X


def test_smacof(filename, target_dimension = 3):
    lengths = np.fromfile(filename, dtype = np.float_)
    l = lengths.shape[0]
    numvertices = int(math.sqrt(l))
    lengths = lengths.reshape((numvertices, numvertices))

    weights = np.vectorize(lambda x: math.exp(-x))(lengths)

    Vplus = find_vplus(weights)

    maximum_iterations = 100
    
    X = smacof_numerical.smacof_embed(target_dimension = target_dimension, weights = weights, lengths = lengths, vplus = Vplus, maximum_iterations = maximum_iterations)
    

    return lengths, weights, Vplus, X
    
