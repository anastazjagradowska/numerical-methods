import math
import numpy as np
import matplotlib.pyplot as plt


def bisection_eigenvalues(A, n, start, end, iteration):
    value = end
    for j in range(iteration):
        count = 0
        w = [0] * (N + 1)
        w[0] = 1
        w[1] = A[0][0] - value
        for i in range(2, len(A) + 1):
            w[i] = (A[i - 1][i - 1] - value) * w[i - 1] - A[i - 1][i - 2]**2 * w[i - 2]
        for i in range(N):
            if w[i] * w[i + 1] < 0:
                count += 1
        if count >= n:
            value = (value + start) / 2
        else:
            value = value + value - start
            start = (value + start) / 2
    return start

def euclidean_norm(A):
    sum = 0
    for i in range(len(A)):
        sum += A[i]*A[i]
    return math.sqrt(sum)

def bisection_eigenvector(A, val):
    x = [0] * len(A)
    x[0] = 1
    x[1] = (val-A[0][0]) / A[0][1]
    for i in range(2, len(A)):
        x[i] = ((val-A[i-1][i-1]) * x[i-1]-A[i-2][i-1] * x[i-2]) / A[i-1][i]
    e = euclidean_norm(x)
    for i in range(len(x)):
        x[i] = x[i] / e
    return x

def gershgorin(A):
    sum = [0] * len(A)
    n = [0] * len(A)
    for i in range(len(A)):
        for j in range(len(A)):
            sum[i] += abs(A[i][j])
        sum[i] -= A[i][i]
        n[i] = A[i][i]
    return -max(sum)-max(n), max(sum)+max(n)

def solutions(A, iter):
    result = []
    vectors = []
    for j in range(1, 6):
        start, end = gershgorin(A)
        x = [0] * len(A)
        y = [0] * len(A)
        for i in range(len(A)):
            x[i] = i + 1
        lam = bisection_eigenvalues(A, j, start, end, iter)
        result.append(lam)
        vectors.append(bisection_eigenvector(A, lam))
    return np.array(result), vectors

def Vpot(x):
    return x**2

a = -5
b = 5
N = 50
x = np.linspace(a, b, 50)
h = x[1] - x[0]

T = np.zeros((N - 2)**2).reshape(N - 2, N - 2)
for i in range(N - 2):
    for j in range(N - 2):
        if i == j:
            T[i, j] = -2
        elif np.abs(i-j) == 1:
            T[i, j] = 1
        else:
            T[i, j] = 0

V = np.zeros((N - 2)**2).reshape(N - 2, N - 2)
for i in range(N - 2):
    for j in range(N - 2):
        if i == j:
            V[i, j] = Vpot(x[i + 1])
        else:
            V[i, j] = 0

H = -T / (2 * h**2) + V

val, vec = solutions(H, N)
vec = np.array(vec)
print(val)
print(vec)
z = np.argsort(val)
z = z[0:5]
energies = (val[z]/val[z][0])
print('Eigenvalues and eigenvectors:')
print(energies)

x_plot = np.linspace(a, b, N)
plt.figure(figsize = (12, 10))

for i in range(len(z)):
    y = []
    y = np.append(y, vec[z[i]])
    y = np.append(y, 0)
    y = np.insert(y, 0, 0)
    plt.plot(x_plot, y, lw=3, label="{} ".format(i))
    plt.xlabel('x', size=14)
    plt.ylabel('$\psi$(x)', size=14)
plt.legend()
plt.title('Wavefunctions', size=14)
plt.show()
