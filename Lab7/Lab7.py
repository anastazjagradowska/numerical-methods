import math

import numpy as np
import matplotlib.pyplot as plt


begin = -5
end = 5

def basic_fun(x):
    return 1/(1+x**2)

def matrix_of_derivative(n,x):
    m=np.zeros((n+1,n+1))

    for i in range(n+1):
        m[i,0]=basic_fun(x[i])

    for j in range(1,n+1):
        for i in range(j,n+1):
            m[i,j]=(m[i,j-1]-m[i-1,j-1])/(x[i]-x[i-j])
    
    return m

def polynomial_interpolation(n,nodes,x):
    der_matrix=matrix_of_derivative(n,nodes)
    acc=0

    for j in range(n+1):
        mult=1
        for i in range(j):
            mult*=(x-nodes[i])
        acc+=der_matrix[j,j]*mult
    
    return acc

def Czebyszew(a,b,n):
    x = [1/2 * ( (a-b)*math.cos(math.pi*((2*i+1)/(2*n+2))) + (a+b)) for i in range(n+1)]
    return np.array(x)

def plot_chart(nodes, title, n):
    x = list(np.linspace(begin, end, num=100))
    y = [polynomial_interpolation(n, nodes, x[j]) for j in range(100)]
    plt.plot(x, y, color='cyan',label="Interpolation")

    y = [basic_fun(x[i]) for i in range(len(x))]
    plt.plot(x, y, color='black',label=r'f(x)=$\frac{1}{1+x^2}$')

    plt.scatter(nodes, [basic_fun(i) for i in nodes], color='red', s=40,label="Interpolation nodes")
    plt.grid()
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend(framealpha=1, frameon=True)
    plt.show()

if __name__ == "__main__":
    nodes=[5,10,15,20]

    for i in nodes:
        for j in range(2):
            if j == 0:
                nodes = np.linspace(-5.0, 5.0, num=i+1)
                title = f"Równoległe n={i}"
            else:
                nodes = Czebyszew(begin, end, i)
                title = f"Czebyszew n={i}"
            
            plot_chart(nodes, title, i)
