import numpy as np
import matplotlib.pyplot as p

print('Podaj N:')
N = int(input())
print('Podaj M:')
M = int(input())
n = N // 2 + M // 2

x_min = -5.0
x_max = 5.0
step = (x_max - x_min) / 100

def func(x):
  return [np.exp(-(i ** 2)) for i in x]

def vector_c(n):
  c = np.zeros(n + 1)
  for k in range(n + 1):
    c[k] = np.math.pow((-1), k) / np.math.factorial(k)
  return c

def vector_A(c):
  A = np.zeros((M // 2, M // 2))
  for i in range(M // 2):
    for j in range(M // 2):
      A[i][j] = c[N // 2 - M // 2 + i + j + 1]
  return A

def vector_y(M, N, c):
  y = np.zeros(M // 2)
  for i in range(M // 2):
    y[i] = -c[N // 2 + 1 + i]
  return y

def calculate_x(A, y):
  x = np.linalg.solve(A, y)
  return x

def calculate_b(M, x):
  b = np.zeros(M // 2 + 1)
  b[0] = 1
  for i in range(M // 2):
    b[M // 2 - i] = x[i]
  return b

def calculate_a(N, c, b):
  a = np.zeros(N // 2 + 1)
  for i in range(N // 2 + 1):
    for j in range(i + 1):
      a[i] += c[i - j] * b[j]
  return a

def vector_R(step, N, a, b):
  R = []
  for i in np.arange(-5, 5 + step, step):
    P = 0
    Q = 0
    for j in range(N // 2 + 1):
      P += a[j] * (i ** (2 * j))
    for j in range(M // 2 + 1):
      Q += b[j] * (i ** (2 * j))
    R.append(P/Q)
  return R

def plot_chart(R, y, step):
  p.plot(np.arange(-5, 5 + step, step), [i for i in R])
  p.plot(np.arange(-5, 5 + step, step), [i for i in y])
  p.legend(["exp($-x^2$)", '$R_{},_{}$'.format(N, M)])
  p.xlabel("x")
  p.ylabel("f(x)")
  p.show()

if __name__ == "__main__":
  print('Rozwiazania:')
  c = vector_c(n)
  A = vector_A(c)
  y = vector_y(M, N, c)
  x = calculate_x(A, y)
  print('x:', x)
  b = calculate_b(M, x)
  print('b:', b)
  a = calculate_a(N, c, b)
  print('a:', a)
  R = vector_R(step, N, a, b)

  x = [i for i in np.arange(-5, 5 + step, step)]
  y = func(x)

  plot_chart(R, y, step)

