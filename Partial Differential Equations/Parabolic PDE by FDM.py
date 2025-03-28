import numpy as np 
import matplotlib.pyplot as plt
import time 
from scipy.sparse import diags
import scipy
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# # Parabolic PDE 

##PDE u_t = u_xx taking h = 1/3 and k = 1/36 
## u(0,t) = 0 = u(1,t) , u(X,0) = sin(pi*x)

# Explicit method
# # Parameters
h = 1/3
k = 1/36
alpha = 1
a = 0
b = 1
M = 15 # Number of spatial steps
N = 2 #time levels 
x = np.linspace(a,b,M+1)
#lambda 
L = (alpha*k)/h**2

#inital condition 
ui = np.sin(np.pi*x) 
#boundary condition 
ui[0] = 0
ui[M] = 0 
# print(ui)
# explcit scheme
u = np.zeros([N+1,M+1])
# print(u)
for j in range(0,N+1):
    for i in range(0,M+1):
        u[0,i] = ui[i]
        u[j,0] = 0
        u[j,M] = 0       

for j in range(0,N):
    for i in range(1,M):          
        u[j+1,i] = u[j,i]+L*(u[j,i-1]-2*u[j,i]+u[j,i+1])
        
print(u)

plt.rcParams['font.family'] = 'Times New Roman'
# Line plot
plt.figure(figsize=(10, 5))
time_levels = [0, 1]  # Time levels to plot (initial and after 1 time step)
for j in time_levels:
    plt.plot(x, u[j], label=f'Time = {j * k:.2f} seconds')

plt.title('Temperature Distribution Over Time')
plt.xlabel('Spatial Domain (x)')
plt.ylabel('Temperature')
plt.legend()
plt.grid()
# plt.savefig('plot 6.png',dpi=500)
plt.show()

# Create a meshgrid for plotting
X, T = np.meshgrid(x, np.linspace(0, N*k, N+1))

# 3D plot
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, T, u, cmap='hot')
ax.set_title('3D Plot of Temperature Distribution')
ax.set_xlabel('Spatial Domain (x)')
ax.set_ylabel('Time (t)')
ax.set_zlabel('Temperature')
# plt.savefig('plot 7.png',dpi=500)
plt.show()



# # _____________________________________________________________________________




# # Implicit method without boundary 

# # Parameters
h = 1/3
k = 1/36
alpha = 1
a = 0
b = 1
m = int((b-a)/h) # Number of spatial steps
n = 2  #time levels 
x = np.linspace(a,b,m+1)
#lambda 
L = (alpha*k)/h**2
# ui = np.sin(np.pi*x) 
# ui = np.array([0.69282032,0.69282032])
p = np.zeros([m-1,m-1])
# print(p)
ai = np.zeros(m-1)
di = np.zeros(m-1)
ui = np.zeros(m-1)
vi = np.zeros(m-1)
# v1= np.zeros(m-1)
# c1 = np.zeros(m-1)
for i in range(0,m-1):
        di[i] = 1+2*L
        ai[i] = -L
        ui[i] = np.sin(np.pi*x[i+1])
print(ui)
# print(ui)
# for z in range():
# for i in range(m-1):
i = 0 
# for i in range(m-1):
while i<n:
    A = diags([ai,di,ai],[-1,0,1],shape = (m-1,m-1)).toarray()
    vi = scipy.sparse.linalg.spsolve(A,ui)
    # vi = np.linalg.solve(A,ui)
    i = i+1
    ui = vi
    # print(A) 
    print(vi)


# #========================================================

# implcit with bundary 
# Parameters
h = 1/3
k = 1/36
alpha = 1
a = 0
b = 1
m = int((b-a)/h) # Number of spatial steps
n = 3  #time levels 
x = np.linspace(a,b,m+1)
# print(x)
#lambda 
L = (alpha*k)/h**2
ui = np.zeros(m+1)
# print(ui)

for i in range(0,m+1):
      ui[i] = np.sin(np.pi*x[i])
      ui[0] = 0
      ui[m] = 0
print(ui)

# ui = np.array([0,0.62712184,0.62712184,0])
p = np.zeros([m-1,m-1])
# print(p)
ai = np.zeros(m+1)
di = np.zeros(m+1)
# b = np.zeros(m+1)
k = 0 
while k<n:
    for i in range(0,m+1):
            di[i] = 1+2*L
            ai[i] = -L
            # ui[i] = np.sin(np.pi*x[i+1])
            # b[0] = 0
            # b[m] = 0
    # for i in range(1,m):
    #     b[i] = ui[i] 
    # print("b",b)
    A = diags([ai,di,ai],[-1,0,1],shape = (m+1,m+1)).toarray()
    # vi = scipy.sparse.linalg.spsolve(A,b)
    vi = np.linalg.solve(A,ui)
    ui = vi
    k = k+1
    # print(A)
    for i in range(0,m+1):
        vi[0] = 0
        vi[m] = 0
    print(vi)


# #__________________________________________________________________________
# Crank Nicolson method

# Parameters
h = 1/3
k = 1/36
alpha = 1
a = 0
b = 1
m = int((b-a)/h) # Number of spatial steps
n = 2#time levels 
x = np.linspace(a,b,m+1)
# print(x)
#lambda 
L = (alpha*k)/h**2
ui = np.zeros(m+1)
# print(ui)

for i in range(0,m+1):
      ui[i] = np.sin(np.pi*x[i])
      ui[0] = 0
      ui[m] = 0
print(ui)

# ui = np.array([0,0.62712184,0.62712184,0])
p = np.zeros([m-1,m-1])
# print(p)
ai = np.zeros(m+1)
di = np.zeros(m+1)
b = np.zeros(m+1)
k = 0 
while k<n:
    for i in range(0,m+1):
            di[i] = 1+L
            ai[i] = -L/2
            # ui[i] = np.sin(np.pi*x[i+1])
            b[0] = 0
            b[m] = 0 
    for i in range(1,m):
            b[i] = L/2*ui[i-1] + (1-L)*ui[i] + L/2*ui[i+1]  
    A = diags([ai,di,ai],[-1,0,1],shape = (m+1,m+1)).toarray()
    # vi = scipy.sparse.linalg.spsolve(A,b)
    vi = np.linalg.solve(A,b)
    ui = vi
    k = k+1
    # print(A)
    for i in range(0,m+1):
        vi[0] = 0
        vi[m] = 0
    print(vi)