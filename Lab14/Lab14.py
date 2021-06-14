import numpy as np
import random
import matplotlib.pyplot as plt

#---1---

x0 = 10
n_1 = 10 ** 4
a = 123
c = 1
m = 2 ** 32

def plot(x_array):
	plt.figure()
	plt.plot(x_array[1::], x_array[0:n_1-1:], "co", markersize=1)

def calculate_1():
	x_array = [x0]
	for i in range(n_1-1):
		x_array.append((a * x_array[i] + c) % m)
	x_array = [x / (m + 1) for x in x_array]

	return x_array

#---2---

u = 4
delta = 3
n = 10**3
k = 10
mi = 4
start = u - delta
end = u + delta
h = (end - start) / k

def dyst(x):
	if x <= mi:
		return((-1 / delta**2) * (-(x ** 2 / 2) + mi * x) + (x / delta) - ((-1 / delta ** 2 ) * (- (mi - delta) ** 2 / 2 + mi * (mi - delta)) + (mi - delta) / delta))
	else:
		return((-1 / delta ** 2) * (x ** 2 / 2 - mi * x) + x / delta - (-1 / delta ** 2 * (mi ** 2 / 2 - mi ** 2) + mi / delta) + 1/2)

def calculate_2():
	Va = 14.0671
	chiSquare_limit = 0
	for _ in range(10000):
		numbers = []

		for i in range(n):
			e1 = random.uniform(0, 1)
			e2 = random.uniform(0, 1)
			x = u + (e1 + e2 - 1) * delta
			numbers.append(x)

		sections = np.zeros(k)
		for number in numbers:
			for j in range(k):
				if number < start + (j + 1) * h:
					sections[j] = sections[j] + 1
					break

		pi_val = [dyst(startSection + h) - dyst(startSection) for startSection in np.arange(start, end, h)]

		chiSquare = sum([(sections[i] - n * pi_val[i])**2 / (n * pi_val[i]) for i in range(k)])

		if chiSquare_limit < Va:
			chiSquare_limit += 1

	return sections, pi_val

def plot_2(sections, pi_val):
	plt.figure()
	plt.bar([i*2 for i in np.arange(start, end, h)], sections/n, align='edge', width=0.8, color="c")
	plt.xticks([i*2-0.2 for i in np.arange(start, end + h, h)], [round(i,1) for i in np.arange(start, end + h, h)])
	plt.plot([i*2+0.4 for i in np.arange(start, end, h)], pi_val, 'r', label='$p_i$')
	plt.show()

if __name__ == "__main__":
	x_array = calculate_1()
	plot(x_array)

	sections, pi_val = calculate_2()
	plot_2(sections, pi_val)

