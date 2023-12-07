#zad2
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange
import numpy as np

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


P_linear = np.zeros(x.shape)
    
for n in range(P):
        if n == 0:
            P_linear += ((data[n+1, 1] - data[n, 1]) / (data[n+1, 0] - data[n, 0]) * (x - data[n, 0])
                             + data[n, 1]) * (x <= data[n+1, 0])
        elif n == P-1:
            P_linear += ((data[n+1, 1] - data[n, 1]) / (data[n+1, 0] - data[n, 0]) * (x - data[n, 0])
                             + data[n, 1]) * (x > data[n, 0])
        else:
            P_linear += ((data[n+1, 1] - data[n, 1]) / (data[n+1, 0] - data[n, 0]) * (x - data[n, 0]) 
                             + data[n, 1]) * ((x > data[n, 0]) & (x <= data[n+1, 0]))

 

lagrange_poly = lagrange(data[:, 0], data[:, 1])

    
plt.plot(data[:,0], data[:,1], 'bo', label = "data points")
plt.plot(x, S, 'b', label = "sklejana stopnia 1")
plt.plot(x, P_linear, 'r', label="sklejana szescienna")
plt.plot(x, lagrange_poly(x), 'g', label="wielomian lagrangea")


plt.xlabel("x")
plt.ylabel("y")
plt.ylim([-4.0, 5.0])
plt.legend(loc=3)

plt.show()