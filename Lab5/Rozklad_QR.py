import copy
import math

import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt

N = 50
A = np.zeros((N, N))
L = -5.0
x = 2 * L / N
iteration = 20

def diagonal_matrix():
    matrix = np.zeros((N, N))
    for i in range(N):
        matrix[i, i] = 1
    return matrix

def euclidean_norm(A):
    sum = 0
    for i in range(N):
        sum += A[0, i] * A[0, i]
    return math.sqrt(sum)

def get_basis(i):
    basis = np.zeros((1, N))
    basis[0, i] = 1
    return basis

def get_H_part(v, i):
    v = np.subtract(v, euclidean_norm(v) * get_basis(i))
    return np.subtract(diagonal_matrix(), (2 / (v @ v.transpose())[0, 0]) * v.transpose() @ v)

def QR_decomposition(A):
    H = []
    v = A.transpose()[0]
    H.append(get_H_part(v, 0))
    R = H[0] @ A

    for i in range(1, N - 1):
        v = copy.copy(R.transpose()[i])
        for j in range(i):
            v[0, j] = 0
        
        H.append(get_H_part(v, i))
        R = H[-1] @ R
    
    Q = diagonal_matrix()
    H.reverse()
    
    for j in H:
        Q = Q @ j

    return R, Q.transpose()

def calculate_eigenvalues(A, P):
    for i in range(iteration):
        R, Q = QR_decomposition(np.asmatrix(A))
        A = R @ Q
        P = P @ Q

    eigenvalues = [A[i, i] for i in range(N)]
    eigenvalues.reverse()

    return P, eigenvalues 

def calculate_eigenvector(N, H, P):
    eigenvector = []
    for i in range(N):
        x = np.full((N, 1), 0)

        for j in range(N - 1, i, -1):
            x[j] = 0
        x[i] = 1

        for j in range(i - 1, -1, -1):
            temp = 0
            for k in range(j + 1, i):
                temp += H[j, k] * x[k]
            x[j] = -temp / (H[j, j] - H[i, i])
        eigenvector.append(P @ x)
    
    eigenvector.reverse()

    return eigenvector

def draw_wavefunctions(eigenvector):
    a = -5
    b = 5
    x_plot = np.linspace(a, b, N)
    plt.figure(figsize=(12, 10))
    for i in range(5):
        y = []
        y = np.append(y, eigenvector[i])
        plt.plot(x_plot, y, lw=3, label="{} ".format(i))
        plt.xlabel('x', size=14)
        plt.ylabel('$\psi$(x)', size=14)
    plt.legend()
    plt.title('Pięć pierwszych funckji falowych', size=14)
    plt.show()

if __name__ == '__main__':
    for i in range(N):
        xi = -L + (i + 1) * x
        A[i, i] = (x ** -2) + (xi ** 2) / 2
        if i != 0:
            A[i, i - 1] = A[i - 1, i] = -1 / (2 * (x ** 2))

    A1 = copy.deepcopy(A)
    P = diagonal_matrix()

    P, eigenvalues = calculate_eigenvalues(A, P)

    print('-- Wartości własne --')
    print(eigenvalues[:5])

    H = inv(P) @ A1 @ P
    eigenvector = calculate_eigenvector(N, H, P)
    draw_wavefunctions(eigenvector)