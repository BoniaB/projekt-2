import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

#zad3
data = np.array(([1.0, 3.0], [2.0, 1.0], [3.5, 4.0], [5.0, 0.0], [6.0, 0.5], [9.0 , -2.0], [9.5, -3.0]))
x = np.linspace(0,9.5, 100)


P = data.shape[0]-1
h = np.zeros(P)
b = np.zeros(P)
u = np.zeros(P)
v = np.zeros(P)
z = np.zeros(data.shape[0]+1)

for i in range(P):
    h[i] = data[i+1,0] - data[i,0]
    b[i] = (6/h[i]) * (data[i+1,0] - data[i,0])

   
for i in range(P):
    u[1] = 2*(h[0] + h[1])
    u[i] = (2*(h[i-1] + h[i])) - ((h[i-1]**2)/u[i-1])
 
    
for i in range(P):    
    v[1] = b[1] - b[0]
    v[i] = (b[i] - b[i-1]) - ((h[i-1]*v[i-1])/u[i-1]) 
    
for i in range(data.shape[0]-1, 0, -1):
    z[i] = (v[i-1] - h[i-1] * z[i+1]) / u[i-1]

z[0] = 0
    

A = np.zeros(data.shape[0])
B = np.zeros(data.shape[0])
C = np.zeros(data.shape[0])

for i in range(data.shape[0]-1):
    A[i] = (z[i+1]-z[i])/(6*h[i])
    B[i] = z[i]/2
    C[i] = (-((h[i]*(z[i+1]+2*z[i]))/6) + (data[i+1,1]-data[i,1])/h[i])


Si = np.zeros(x.shape)
for i in range(data.shape[0]-1):
    Si += (data[i,1] + (x-data[i,0])*(C[i] + (x-data[i,0])*(B[i] + ((x-data[i,0])*A[i])))) * ((x > data[i,0]) & (x <= data[i+1,0]))


S = np.zeros(x.shape)

for i in range(data.shape[0]-1):
    if i == 0:
        S += Si * (x <= data[i+1, 0])
    elif i == data.shape[0]-2:
        S += Si * (x > data[i, 0])
    else:
        S += Si * ((x > data[i, 0]) & (x <= data[i+1, 0]))
    



ti = np.array([1.0, 2.0, 3.5, 5.0, 6.0, 9.0, 9.5])
yi = np.array([3.0, 1.0, 4.0, 0.0, 0.5, -2.0, -3.0])


cs = CubicSpline(data[:,0], data[:,1], bc_type='natural')

plt.plot(data[:,0], data[:,1], 'bo', label = "data points")
plt.plot(x, S, 'r', label = "sklejana szescienna")
plt.plot(x, cs(x), 'g', label='Cubic Spline')

plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc=1)

plt.show()