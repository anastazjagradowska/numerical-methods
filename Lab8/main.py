import matplotlib.pyplot as plt
import numpy as np
import math
from numpy.linalg import inv

n = [5, 8, 15, 21]
a = -5
b = 5


def fun1(x):
    return 1 / (1 + x ** 2)


def fun2(x):
    return math.cos(2 * x)


def derivative_second(x, fun, dx=0.01):
    return (fun(x - dx) - 2 * fun(x) + fun(x + dx)) / dx ** 2


def d_M(nodes, fun):
    h = []
    for i in range(1, len(nodes)):
        h.append(nodes[i] - nodes[i - 1])
    lambd = []
    mi = []
    for i in range(len(h) - 1):
        lambd.append(h[i + 1] / (h[i] + h[i + 1]))
        mi.append(1 - h[i + 1] / (h[i] + h[i + 1]))

    d = [0] * len(nodes)
    for i in range(1, len(nodes) - 1):
        d[i] = (6 / (h[i - 1] + h[i]) * (
                (fun(nodes[i + 1]) - fun(nodes[i])) / h[i] - (fun(nodes[i]) - fun(nodes[i - 1])) / h[i - 1]))

    A = np.eye(len(nodes))
    for i in range(len(nodes)):
        A[i, i] += 1
    for i in range(len(nodes) - 2):
        A[i + 1, i + 2] = lambd[i]
        A[i + 1, i] = mi[i]

    d1 = [0] * len(nodes)
    for i in range(1, len(nodes) - 1):
        d1[i] = 3
    m = d @ np.linalg.inv(A)
    return m


def d_Sx(nodes, fun, x):
    h = []
    for k in range(1, len(nodes)):
        h.append(nodes[k] - nodes[k - 1])

    i = 0
    for k in range(len(nodes) - 1):
        i = k + 1 if nodes[k] < x < nodes[k + 1] else i
        i = 1 if x <= min(nodes) else i
        i = len(nodes) - 1 if x >= max(nodes) else i

    m = d_M(nodes, fun)
    A = (fun(nodes[i]) - fun(nodes[i - 1])) / h[i - 1] - h[i - 1] / 6 * (m[i] - m[i - 1])
    B = fun(nodes[i - 1]) - m[i - 1] * h[i - 1] ** 2 / 6

    return m[i - 1] * (nodes[i] - x) ** 3 / (6 * h[i - 1]) + m[i] * (x - nodes[i - 1]) ** 3 / (6 * h[i - 1]) + A * (
            x - nodes[i - 1]) + B


def calculation_error(y1, y2):
    e = 0
    for i in range(len(y1)):
        e = abs(y1[i] - y2[i]) if (abs(y1[i] - y2[i]) > e) else e
    return e


def plot_interpolation(j, fun, label):
    nodes = np.linspace(a, b, num=j)
    x = list(np.linspace(a, b, num=100))
    y1 = [d_Sx(nodes, fun, i) for i in x]
    plt.plot(x, y1, color='cyan', label="Interpolation")
    plt.scatter(nodes, [fun(i) for i in nodes], color='hotpink', s=40, label="Interpolation nodes")
    y2 = [fun(x[i]) for i in range(len(x))]
    plt.plot(x, y2, color='black', label=label)
    plt.grid()
    plt.title(
        f"Equidistant nodes n={j}, h={round(nodes[1] - nodes[0], 2)}  \u03B5\u2098\u2090\u2093={round(calculation_error(y1, y2), 2)}")
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend(framealpha=1, frameon=True)
    plt.show()


def plot_derivative_second():
    nodes = np.linspace(a, b, 11)
    m = d_M(nodes, fun1)
    x = list(np.linspace(a, b, num=100))
    derivative = [derivative_second(nodes[i], fun1) for i in range(len(m))]
    y1 = [d_Sx(nodes, fun2, i) for i in x]
    plt.plot(nodes, derivative, color='black')
    plt.scatter(nodes, m, color='cyan', s=40, label="m derivative")
    plt.scatter(nodes, derivative, color='hotpink', s=40, label="approx derivative")
    y2 = [fun2(x[i]) for i in range(len(x))]
    plt.grid()
    plt.title(
        f"Second derivatives n={11}, h={round(nodes[1] - nodes[0], 2)}, \u03B5\u2098\u2090\u2093={round(calculation_error(y1, y2), 2)}")
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend(framealpha=1, frameon=True)
    plt.show()


if __name__ == "__main__":
    for j in n:
        plot_interpolation(j, fun1, r'f(x)=$\frac{1}{1+x^2}$')
        plot_interpolation(j, fun2, r'f(x)=cos(2x)$')
    plot_derivative_second()
