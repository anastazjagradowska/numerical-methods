import numpy as np
import math
from matplotlib import pyplot as pl
import time


def euclides_norm(array):
    sum = 0
    for i in range(len(array)):
        sum += (array[i] ** 2)
    return math.sqrt(sum)


def method_of_steepest_descent(A, b):
    N = len(A)
    iteration_list = []
    x_norm_list = []
    r_norm_list = []
    alpha_list = []
    x = np.zeros((N, 1))
    i = 1

    while True:
        r = b - A @ x
        alpha = (np.matmul(r.transpose(), r))[0, 0] / (np.matmul(np.matmul(r.transpose(), A), r))[0, 0]
        alpha_list.append(alpha)
        x = x + alpha * r
        iteration_list.append(i)
        x_norm_list.append(euclides_norm(x))
        r_norm_list.append(euclides_norm(r))
        i += 1
        if euclides_norm(r) < 1e-3:
            break

    return x, iteration_list, x_norm_list, r_norm_list, alpha_list


def method_of_steepest_descent_to_file_save(iteration_list, r_norm_list, alpha_list, x_norm_list):
    file = open("method_of_steepest_descent" + ".txt", "w")
    for i in range(len(iteration_list)):
        file.write(str(iteration_list[i]) + "\t" + str(r_norm_list[i]) + "\t" + str(alpha_list[i]) + "\t" +
                   str(x_norm_list[i]) + "\n")
    file.close()


def method_of_conjugat_to_file_save(iteration_list, r_norm_list, alpha_list, beta_list, x_norm_list):
    file = open("method_of_conjugat_descent" + ".txt", "w")
    for i in range(len(iteration_list)):
        file.write(str(iteration_list[i]) + "\t" + str(r_norm_list[i]) + "\t" + str(alpha_list[i]) + "\t" +
                   str(beta_list[i]) + "\t" + str(x_norm_list[i]) + "\n")
    file.close()


def show_method_of_conjugate_diagram(iteration_list, x_norm_list, r_norm_list):
    show_diagram(iteration_list, x_norm_list, r_norm_list,
                 "Wartości norm w kolejnych iteracjach metodą sprężonego gradientu")


def show_method_of_steepest_descent_diagram(iteration_list, x_norm_list, r_norm_list):
    show_diagram(iteration_list, x_norm_list, r_norm_list,
                 "Wartości norm w kolejnych iteracjach metodą największego spadku")


def show_diagram(iteration_list, x_norm_list, r_norm_list, diagram_title):
    F_x, ax1 = pl.subplots()

    pl.title(diagram_title)

    ax1.set_xlabel('numer iteracji')
    ax1.set_ylabel('norma wektora rozwiązań', color='darkred')
    ax1.plot(iteration_list, x_norm_list, 'r.--', label='norma ||x||', color="pink")
    ax1.tick_params(axis='y')

    ax2 = ax1.twinx()

    ax2.plot(iteration_list, r_norm_list, 'g.--', label='norma ||r||', color="k")
    ax2.tick_params(axis='y')
    ax2.set_ylabel('norma wektora reszt', color='deepskyblue')
    ax2.set_yscale("log")

    F_x.tight_layout()
    pl.legend()
    pl.show()


def conjugate_method(A, b):
    n = len(b)
    x = np.ones(n)
    r = np.dot(A, x) - b
    p = - r
    r_k_norm = np.dot(r, r)
    i = 0
    x_norm_list = []
    r_norm_list = []
    iteration_list = []
    alpha_list = []
    beta_list = []
    while euclides_norm(r) > ((10 ** (-6))):
        i = i + 1
        Ap = np.dot(A, p)
        alpha = r_k_norm / np.dot(p, Ap)
        alpha_list.append(alpha)
        x += alpha * p
        r += alpha * Ap
        x_norm_list.append(euclides_norm(x))
        r_norm_list.append(euclides_norm(r))
        iteration_list.append(i)
        r_kplus1_norm = np.dot(r, r)
        beta = r_kplus1_norm / r_k_norm
        beta_list.append(beta)
        r_k_norm = r_kplus1_norm
        p = beta * p - r
    return x, iteration_list, x_norm_list, r_norm_list, alpha_list, beta_list


if __name__ == '__main__':

    N = 1000
    M = 5.0
    eps = 1e-3
    A = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if np.abs(i - j) <= M:
                A[i, j] = (1 / (1 + np.abs(i - j)))

    b = np.arange(0.0, N).reshape(N, 1)

    start_time = time.time()

    x, iteration_list, x_norm_list, r_norm_list, alpha_list = method_of_steepest_descent(A, b)

    show_method_of_steepest_descent_diagram(iteration_list, x_norm_list, r_norm_list)

    print("Metoda największego spadku: " + str((time.time() - start_time)) + " sekund")

    method_of_steepest_descent_to_file_save(iteration_list, r_norm_list, alpha_list, x_norm_list)

    start_time = time.time()
    b = []
    for i in range(N):
        b.append(i)

    x, iteration_list, x_norm_list, r_norm_list, alpha_list, beta_list = conjugate_method(np.array(A), np.array(b))
    method_of_conjugat_to_file_save(iteration_list, r_norm_list, alpha_list, beta_list, x_norm_list)
    print("Metodą sprężonego gradientu: " + str((time.time() - start_time)) + " sekund")
    show_method_of_conjugate_diagram(iteration_list, x_norm_list, r_norm_list)
