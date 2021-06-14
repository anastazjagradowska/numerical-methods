import matplotlib.pyplot as plt
import numpy as np
import math

def f(x):
	return x / (4* x**2 + 1)

def fc2(x):
	return x**10

def fc3x(x):
	return np.sin(x) ** 2

def fc3y(x):
	return np.sin(x) ** 4

def check(x):
	return (1/(2*2**2)) * math.log( 2**2 * x**2 + 1**2 )

def legendre():
	IntC1 = []
	for n in range(2,21):
		x, w = np.polynomial.legendre.leggauss(n)
		s = sum(w*f(x+1))
		print(s)
		IntC1.append(s)
	IntC1 = np.array(IntC1)
	return IntC1

def gauss_Laguerre_quadrature(fc2):
	IntC2 = []
	for n in range(2,21):
		x,w = np.polynomial.laguerre.laggauss(n)
		s = sum(w*fc2(x))
		IntC2.append(s)
	IntC2 = np.array(IntC2)
	return IntC2

def hermgauss(fc3x, fc3y):
	IntHer = []
	for n in range(2,16):
		x, w = np.polynomial.hermite.hermgauss(n)
		sx = sum(w * fc3x(x))
		sy = sum(w * fc3y(x))
		IntHer.append(sx * sy)
	IntHer = np.array(IntHer)
	return IntHer

def plot(Int_tab, fc, color, range_max):
	plt.figure()
	plt.plot([i for i in range(2, range_max)], np.abs(Int_tab - fc), color)
	plt.yscale('log')
	plt.grid()


if __name__ == '__main__':
	IntC1 = legendre()
	IntC2 = gauss_Laguerre_quadrature(fc2)
	IntHer = hermgauss(fc3x, fc3y)
	plot(IntC1, check(2) - check(0), '-or', 21)
	plot(IntC2, math.factorial(10), '-og', 21)
	plot(IntHer, 0.1919832644, '-ob', 16)
	plt.show()
