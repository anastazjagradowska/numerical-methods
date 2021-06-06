import numpy as np
import math
import matplotlib.pyplot as plt

a = 0
b = math.pi
p_array = [5, 10, 25, 50, 100]
i = 30
# m = 0
# k = 1
# m = 1
# k = 1
m = 5
k = 5

def fun(x, m, k):
    return x ** m * math.sin(k * x)

def sum(i, m, k, x):
    return (-1) ** i * ((k * x) ** (2 * i + m + 2)) / (
            k ** (m + 1) * math.factorial(2 * i + 1) * (2 * i + m + 2))

def getSum(numberOfIterations, m, k, a, b):
    sums = []
    I = 0
    for i in range(numberOfIterations):
        I = I + sum(i, m, k, b) - sum(i, m, k, a)
        sums.append(I)
    return sums

def Simpson(p_array, m, k, a, b):
    I_array = []
    for p in p_array:
        h = (b - a) / p
        I = 0
        for i in np.arange(a, b, h):
            left = i
            middle = i + h / 2
            right = i + h
            Ii = 1 / 3 * h / 2 * (1 * fun(left, m, k) + 4 * fun(middle, m, k) + 1 * fun(right, m, k))
            I = I + Ii
        I_array.append(I)
    return I_array

def plot_sum_value_changes(sum_array):
    plt.plot(range(1, i + 1), sum_array)
    plt.plot(range(1, i + 1), sum_array, "bo")
    plt.xlabel('$i$')
    plt.ylabel('$I$')
    plt.title('k = {}, m = {}'.format(k, m))
    plt.show()

def plot_dependence_on_the_number_of_nodes(simpson):
    plt.plot([2 * p + 1 for p in p_array], [abs(C-I) for C in simpson])
    plt.plot([2 * p + 1 for p in p_array], [abs(C-I) for C in simpson], "bo")
    plt.xlabel('$n$')
    plt.ylabel('$|C-I|$')
    plt.yscale("log")
    plt.title('k = {}, m = {}'.format(k, m))
    plt.show()

if __name__ == "__main__":
    sum_array = getSum(i, m, k, a, b)
    I = sum_array[-1]
    simpson = Simpson(p_array, m, k, a, b)
    plot_sum_value_changes(sum_array)
    plot_dependence_on_the_number_of_nodes(simpson)
