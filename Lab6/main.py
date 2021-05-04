import numpy as np
import matplotlib.pyplot as p

eps = 1e-6

def get_delta(x, y):
    M = np.array([[(2 * y ** 2 - 6 * x * y), (4 * x * y - 3 * x ** 2)],
                  [(2 * x * (y ** 3) + 2 * y), (3 * (x ** 2) * (y ** 2) + 2 * x)]])
    v = np.array([[2 * x * y ** 2 - 3 * (x ** 2) * y - 2], [(x ** 2) * (y ** 3) + 2 * x * y - 12]])
    M = np.linalg.inv(M)
    v = -M @ v
    return v.transpose()

def calculate_vectors(r):
    delta_X = []
    iteration = 0
    while 1:
        print(iteration, r[-1])
        d = get_delta(r[-1][0][0], r[-1][0][1])
        r.append(r[-1] + d)
        delta_X.append(np.linalg.norm(d))
        iteration += 1
        if delta_X[-1] < eps:
            break
    return delta_X

def plot_vectors(r):

    p.plot(calculate_vectors(r), marker='o', markersize=3)
    p.xlabel('Numer iteracji')
    p.ylabel('\u0394 r')
    p.yscale('log')
    p.grid()

    p.figure()
    for i in range(len(r) - 1):
        dx = r[i + 1][0][0] - r[i][0][0]
        dy = r[i + 1][0][1] - r[i][0][1]
        p.quiver(r[i][0][0], r[i][0][1], dx, dy, angles='xy', scale_units='xy', scale=1)
    p.xlim(0, 12)
    p.ylim(-6, 20)
    p.plot(r[-1][0][0], r[-1][0][1], marker='o', markersize=3, color='cyan')
    p.xlabel('x')
    p.ylabel('y')
    p.grid()
    p.show()

r = [np.array([[10, 10]])]
plot_vectors(r)

r = [np.array([[10, -4]])]
plot_vectors(r)
