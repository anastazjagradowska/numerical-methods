import numpy as np
import math
import random
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt

k = 12
Nk = 2 ** k
T = 1.0
t_max = 3 * T
dt = t_max / Nk
t = np.arange(0, t_max, dt)
sigma = T / 20
omega = 2 * math.pi / T

def F0(t):
  return math.sin(1 * omega * t) + math.sin(2 * omega * t) + math.sin(3 * omega * t)

def F(t):
  return F0(t) + random.uniform(-1/2, 1/2)

def G(t):
  return 1 / (sigma * math.sqrt(2 * math.pi)) * math.exp(- t**2 / (2 * sigma**2))

def calculate_Fourier():
    f0_array = np.asarray([F0(t) for t in t])
    f_array = np.asarray([F(t) for t in t])
    g1_array = np.asarray([G(t) for t in t])
    g2_array = np.asarray([G(-t) for t in t])

    g_array = g2_array + g1_array

    fFourier = fft(f_array)
    gFourier = fft(g_array)
    fgFourier = np.multiply(fFourier, gFourier)

    fgIFourier = ifft(fgFourier)
    f_max = np.amax(np.asarray([abs(value) for value in fgIFourier]))

    return f_array, f0_array, fgIFourier, f_max

def plot(f_array, f0_array, fgIFourier, f_max):
    plt.plot(t, f_array, label='$f$', color="cyan")
    plt.plot(t, f0_array, label='$f_0$', color="fuchsia")
    plt.plot(t, fgIFourier * 2.5 / f_max, label='$wygladzone$', color="blue")
    plt.xlabel('$t$')
    plt.ylabel('$f(t)$')
    plt.title('k = {}'.format(k))
    plt.legend()
    plt.show()

if __name__ == "__main__":
    f_array, f0_array, fgIFourier, f_max = calculate_Fourier()
    plot(f_array, f0_array, fgIFourier, f_max)
