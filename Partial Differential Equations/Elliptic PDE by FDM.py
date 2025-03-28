import numpy as np 
import matplotlib.pyplot as plt
import time 
from scipy.sparse import diags
import scipy
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

#Elliptic PDE 

# 1. Laplace Equation u_xx + u_yy = 0 
# # # # u(x, 0) = 10, u(0, y) = 10, u(x, 1) = 20, u(1, y) = 20

h = 0.25
k = 0.25
x0 = 0
xn = 1
y0 = 0
yn = 1
M = int((xn-x0)/h)
N = int((yn-y0)/k)
x = np.linspace(0,1,M+1)
y = np.linspace(0,1,N+1)
u = np.zeros([M+1,N+1])
# print(u)
# boundary condition 
for k in range(20):
    for i in range(0,M+1):
        for j in range(0,N+1):
            u[0,j] = 10
            u[i,0] = 10
            u[N,j] = 20
            u[i,M] = 20
    # print(u)
    for i in range(1,M):
        for j in range(1,N):
            u[i,j] = 1/4*(u[i+1,j]+u[i-1,j]+u[i,j+1]+u[i,j-1])
print("1",u)

# Plotting the solution
X, Y = np.meshgrid(x, y)
plt.figure(figsize=(8, 6))
contour = plt.contourf(X, Y, u.T, levels=50, cmap='viridis',)
plt.colorbar(contour)
plt.title('Solution of Laplace Equation (u_xx + u_yy = 0)')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.show()

# Plotting the solution as a linear graph
plt.figure(figsize=(8, 6))
y_fixed = int(N / 2)  # Choose a fixed y value (middle row)
plt.plot(x, u[:, y_fixed], marker='o', label=f'u(x, y={y[y_fixed]})')
plt.title('Linear Graph of the Laplace Equation Solution')
plt.xlabel('x')
plt.ylabel('u(x, y)')
plt.grid()
plt.legend()
plt.show()

# #______________________________________________________________________________

# # 2. Poission's Equation u_xx + u_yy = sin(x.y)

h = 0.25
k = 0.25
x0 = 0
xn = 1
y0 = 0
yn = 1
M = int((xn-x0)/h)
N = int((yn-y0)/k)
x = np.linspace(0,1,M+1)
y = np.linspace(0,1,N+1)
u = np.zeros([M+1,N+1])
# boundary condition 
for k in range(30):
    for i in range(0,M+1):
        for j in range(0,N+1):
            u[0,j] = 5
            u[i,0] = 10
            u[N,j] = 15
            u[i,M] = 20
            # print(u)
    for i in range(1,M):
        for j in range(1,N):
            u[i,j] = 1/4*(u[i+1,j]+u[i-1,j]+u[i,j+1]+u[i,j-1])-h**2/4*(np.sin(x[i]*y[j]))
print("2",u)

# Plotting the solution
X, Y = np.meshgrid(x, y)
plt.figure(figsize=(8, 6))
contour = plt.contourf(X, Y, u.T, levels=50, cmap='viridis',)
plt.colorbar(contour)
plt.title('Solution of Laplace Equation (u_xx + u_yy = x^2+y^2)')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.show()

# Plotting the solution as a linear graph
plt.figure(figsize=(8, 6))
y_fixed = int(N / 2)  # Choose a fixed y value (middle row)
plt.plot(x, u[:, y_fixed], marker='o', label=f'u(x, y={y[y_fixed]})')
plt.title('Solution of Laplace Equation (u_xx + u_yy = x^2+y^2)')
plt.xlabel('x')
plt.ylabel('u(x, y)')
plt.grid()
plt.legend()
plt.show()





# #______________________________________________________________________________

# # Solution of Partial Dierential Equation
# # Elliptic PDE 
# # Laplace Equation u_xx + u_yy = 0 
# # F(X,Y)=x**2y**2

h = 1
k = 1
x0 = 0
xn = 4
y0 = 0
yn = 4
M = int((xn-x0)/h)
N = int((yn-y0)/k)
x = np.linspace(0,1,M+1)
y = np.linspace(0,1,N+1)
u = np.zeros([M+1,N+1])
ti = time.time()
# boundary condition 
for k in range(10):
    for i in range(0,M+1):
        for j in range(0,N+1):
            u[0,j] = 0
            u[i,0] = 0
            u[N,1] = 16
            u[N,2] = 64
            u[N,3] = 144
            u[N,N] = 256
            u[1,M] = 16
            u[2,M] = 64
            u[3,M] = 144
# print(u)
    for i in range(1,M):
        for j in range(1,N):
            u[i,j] = 1/4*(u[i+1,j]+u[i-1,j]+u[i,j+1]+u[i,j-1])
            print(u)
tf = time.time()
T = tf-ti
print(T)



# _______________________________________________________________________________

iter = 10
a=0
b=4
h=1
k=1
nx=int((b-a)/h)
ny=int((b-a)/k)

x=np.linspace(a,b,nx+1)
y=np.linspace(a,b,ny+1)
 
u = np.zeros([ny+1,nx+1])

def f(x,y):
    return (x*y)**2
ti = time.time()
for _ in range(iter):
    for i in range(nx+1):
        for j in range(ny+1):
            u[0,j]=0
            u[nx,j]=f(x[nx],y[j])
            u[i,0]=0
            u[i,ny]=f(x[i],y[ny])
    for i in range(1,nx):
        for j in range(1,ny):
            u[i,j]= (u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1])/4
tf = time.time()
T = tf-ti
print(T)

# Plot the results
X, Y = np.meshgrid(x, y)
plt.contourf(X, Y, u, cmap='viridis')
plt.colorbar(label="Solution Value (u)")
plt.title("Contour Plot of u(x, y)")
plt.xlabel("x")
plt.ylabel("y")
plt.show()

