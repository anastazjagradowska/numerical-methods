import numpy as np
import matplotlib.pyplot as p

x_min = -1.5
x_max = 1

h = 0.01
# n = 10
n = 100
eps = 1e-8

# x1 = -0.5
# x1 = -0.9
x1 = 1.5
x2 = x1 + h
x3 = x2 + h


def func1(x):
    return [np.log(pow(i, 5) + 3 * pow(i, 2) + i + 9) for i in x]


def func2(x):
    return [pow(i, 6) for i in x]


def iterable_function(x):
    return np.log(pow(x, 5) + 3 * pow(x, 2) + x + 9)


def iterable_function2(x):
    return pow(x, 6)


def F2(x1, x2, f):
    return (f(x2) - f(x1)) / (x2 - x1)


def F3(x1, x2, x3, f):
    return (F2(x2, x3, f) - F2(x1, x2, f)) / (x3 - x1)


def approximate_position_of_xm(x1, x2, x3, f):
    return (x1 + x2) / 2 - (F2(x1, x2, f) / (2 * F3(x1, x2, x3, f)))


def powell(approximate_position_of_xm, iter_func, F2, eps, x1, x2, x3):
    xm = []
    F_2 = []
    F_3 = []
    array_x1 = [x1]
    array_x2 = [x2]
    array_x3 = [x3]
    for i in range(n):
        # f = iterable_function
        f = iterable_function2
        m = approximate_position_of_xm(x1, x2, x3, f)
        xm.append(m)
        F_2.append(F2(x1, x2, f))
        F_3.append(F3(x1, x2, x3, f))

        if (np.fabs(x1 - m) < eps):
            break
        if (np.fabs(x2 - m) < eps):
            break
        if (np.fabs(x3 - m) < eps):
            break

        if (np.fabs(m - x1) > np.fabs(m - x2) and np.fabs(m - x1) > np.fabs(m - x3)):
            x1 = m
        elif (np.fabs(m - x2) > np.fabs(m - x1) and np.fabs(m - x2) > np.fabs(m - x3)):
            x2 = m
        else:
            x3 = m

        array_x1.append(x1)
        array_x2.append(x2)
        array_x3.append(x3)
        print("x1 = ", x1, "  x2 = ", x2, "  x3 = ", x3, "  x_m = ", m)

        print(array_x1, array_x2, array_x3)
    return xm, F_2, F_3, array_x1, array_x2, array_x3


def draw_chart(func, iter_func, x_min, x_max, h):
    analytical_values_of_x = [i for i in np.arange(x_min, x_max + h, h)]
    analytical_values_of_y = func(analytical_values_of_x)

    p.plot([i for i in analytical_values_of_x], [i for i in analytical_values_of_y], linestyle = '-', marker = '',
           color = 'gainsboro')
    p.xlabel("x")
    p.ylabel("f(x)")
    p.plot(xm[len(xm) - 1], iter_func(xm[len(xm) - 1]), marker = 'o', color = 'mediumturquoise')
    p.plot(xm[0], iterable_function(xm[0]), marker = 'o', color = 'cadetblue')
    p.show()


def draw_difference_quotient(F_2, F_3):
    p.plot([i for i in range(len(F_2))], [i for i in F_2], linestyle = '-', marker = '.')
    p.plot([i for i in range(len(F_3))], [i for i in F_3], linestyle = '-', marker = '.')
    p.xlabel('k')
    p.ylabel('Wartość ilorazu różnicowego')
    p.show()


def draw_another_approximation(array_x1, array_x2, array_x3, xm):
    p.plot([i for i in range(len(array_x1))], [i for i in array_x1], linestyle = '-', marker = '.')
    p.plot([i for i in range(len(array_x2))], [i for i in array_x2], linestyle = '-', marker = '.')
    p.plot([i for i in range(len(array_x3))], [i for i in array_x3], linestyle = '-', marker = '.')
    p.plot([i for i in range(len(xm))], [i for i in xm], linestyle = '-', marker = '.')
    p.xlabel('k')
    p.ylabel('Wartość kolejnego przybliżenia')
    p.show()


if __name__ == '__main__':
    xm, F_2, F_3, array_x1, array_x2, array_x3 = powell(approximate_position_of_xm, iterable_function2, F2, eps, x1, x2,
                                                        x3)

    draw_chart(func2, iterable_function2, x_min, x_max, h)
    draw_difference_quotient(F_2, F_3)
    draw_another_approximation(array_x1, array_x2, array_x3, xm)
