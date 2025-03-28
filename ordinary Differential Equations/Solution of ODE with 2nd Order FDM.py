import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
import math

# # # BVP y'' = sinx with y(0) = y(1) = 0
# # given parameters
a=0 #x_0
b=1 #xn
n= 100 #number of points 
p= 0 #Coieefficent of y'
q = 0 #Coieefficent of y
h=(b-a)/n #step Size
print(h)
# values of y 
alpha  = 0 #y at x = 0 
beta = 0 #y at x = n 
# print(h)
#x values
x=np.linspace(a,b,n+1)
# Defintion of array for diagonals 
li = np.zeros(n-1) #lower 
di = np.zeros(n-1) # middel
ui = np.zeros(n-1) #upper

#Right hand side of system of equation
gamma=np.zeros(n-1)

#setting Boundary condition in Right hand side
for i in range(1,n-2):
    gamma[i]= np.sin(x[i+1])
gamma[0] = np.sin(x[1])-li[0]*alpha
gamma[n-2] = np.sin(x[n-1])-li[0]*alpha
# print("rhs",gamma)

# Formation of diagonal
for i in range(2, n):
    li[i-2] = 1/h**2-(p/2*h)
for i in range(1, n):
    ui[i-1] = 1/h**2+(p/2*h)
for i in range(0, n-1):
    di[i] = -2/(h**2)+q
    
# Sparse matrix A setup
A = diags([li, di, ui], [-1,0,1], shape=(n-1, n-1)).toarray()
# print(A)

# Solve the linear system
# rslt = spsolve(A,gamma)
rslt = np.linalg.solve(A,gamma)
# print('rslt=',rslt)


# Exact solution
u_exact= np.zeros(n-1)
for i in range(1,n):
    # u_exact[i-1] = x[i]*np.sin(x[i])-np.sin(x[i]) 
    u_exact[i-1] = (1/2)*(1/np.tan(1))*np.sin(x[i])-(1/2)*x[i]*np.cos(x[i])

# print('u_exact=',u_exact)

# Calculate and print the errors
errors = max(np.abs(rslt - u_exact))
print(errors)

# for i in range(n-1):
    # print(f"x = {x[i+1]:.2f}, Numerical = {rslt[i]:.4f}, Exact = {u_exact[i]:.4f}, Error = {errors[i]:.4e}")


# Plot the numerical and exact solutions
plt.plot(x[:n-1], rslt, '-o', label='Numerical Solution (FD)')
plt.plot(x[:n-1], u_exact, '--', label='Exact Solution')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.title('Comparison of Numerical and Exact Solutions')
plt.legend()
plt.grid(True)
# plt.savefig('Q1.png',dpi=500)
plt.show()

#___________________________________________________________________________________________________________________________

# BVP y''+ xy'-y = xcosx-2sinx  with y(0) =0,  y(1) = sin1

import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
import math

# given parameters
a=0
b= 1
n= 4
h=(b-a)/n
print(h)

# boundary value conditions
x=np.linspace(a,b,n+1)
# print(x)
# tridiagonal sys
li=np.zeros(n-2)
di=np.zeros(n-1)
ui=np.zeros(n-1)


for i in range(2, n):
    li[i-2] = (1/h**2)-x[i]/(2*h)
    # li[i-2] = (1)-x[i]*h/(2)
for i in range(1, n):
    ui[i-1] = (1/h**2)+x[i]/(2*h)    
    # ui[i-1] = (1)+x[i]*h/(2)
    
for i in range(0, n-1):
    # li[i] = (1/h**2)-x[i]/(2*h)
    di[i] = -2/(h**2)-1
    # di[i] = -2-h**2
    
# Sparse matrix A setup
A = diags([li, di, ui], [-1,0,1], shape=(n-1, n-1)).toarray()
# print(A)

# Solve the linear system
# rslt = spsolve(A,gamma)
gamma=np.zeros(n-1)

for i in range(1,n-1):
    gamma[i-1]= x[i]*np.cos(x[i])-2*np.sin(x[i])
    gamma[n-2] =  (x[n-2]*np.cos(x[n-2])-2*np.sin(x[n-2]))-ui[n-2]*np.sin(1) 
# print(gamma)

rslt = np.linalg.solve(A,gamma)
# print('rslt=',rslt)



# Exact solution
u_exact= np.zeros(n-1)
for i in range(1,n):
    # u_exact[i] = x[i]*(1-x[i])/2
    u_exact[i-1] = np.sin(x[i])

# print('u_exact=',u_exact)

# Calculate and print the errors
errors = max(np.abs(rslt - u_exact))
print(errors)

# for i in range(n-1):
    # print(f"x = {x[i+1]:.2f}, Numerical = {rslt[i]:.4f}, Exact = {u_exact[i]:.4f}, Error = {errors[i]:.4e}")


# Plot the numerical and exact solutions
plt.plot(x[:n-1], rslt, '-o', label='Numerical Solution (FD)')
plt.plot(x[:n-1], u_exact, '--', label='Exact Solution')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.title('Comparison of Numerical and Exact Solutions')
plt.legend()
plt.grid(True)
plt.savefig('Q2.png',dpi=500)
plt.show()

#__________________________________________________________________________________________________________

# y'' = -1 y(0)=0=y(1)

# import numpy as np
# from scipy.sparse import diags
# import matplotlib.pyplot as plt
# from scipy.sparse.linalg import spsolve
# import math

# given parameters
a=0
b=1
n=4
h=(b-a)/n
print(h)
p=0
q=0
# boundary value conditions
x=np.linspace(a,b,n+1)
print(x)
# tridiagonal sys
li=np.zeros(n-2)
di=np.zeros(n-1)
ui=np.zeros(n-2)
gamma=np.zeros(n-1)

for i in range(0,n-1):
  gamma[i]=-h**2
# print(x[i+1])
print(gamma)
# gamma = np.exp(x)

# for i in range(0, n-1): wrong 
li[:] = 1-h*p/2
di[:] = -2+(h**2)+q
ui[:] = 1+h*p/2

# Sparse matrix A setup
A = diags([li, di, ui], [-1,0,1], shape=(n-1, n-1)).toarray()
print(A)

# Solve the linear system
rslt = spsolve(A,gamma)
print('rslt=',rslt)


# Exact solution
u_exact= np.zeros(n-1)
for i in range(n-1):
    u_exact[i] = x[i]*(1-x[i])/2

print('u_exact=',u_exact)

# Calculate and print the errors
errors = np.abs(rslt - u_exact)
for i in range(n-1):
    print(f"x = {x[i]:.2f}, Numerical = {rslt[i]:.4f}, Exact = {u_exact[i]:.4f}, Error = {errors[i]:.4e}")


# Plot the numerical and exact solutions
plt.plot(x[:n-1], rslt, '-o', label='Numerical Solution (FD)')
plt.plot(x[:n-1], u_exact, '--', label='Exact Solution')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.title('Comparison of Numerical and Exact Solutions')
plt.legend()
plt.grid(True)
plt.show()

#__________________________________________________________________________________________________________________

# # linear ODE constant Coeficient 
# # left dirichlet and right neumann BVP
# # y'' = e^x with y(0)=1 and y'(1)=1 Exact Solution e^x(1-e)x
# # given parameters
a=0
b=1
yL = 1
yR =1 
n= 3
h=(b-a)/n
print(h)
p=0
q=0
# boundary value conditions
x=np.linspace(a,b,n+1)
print(x)
# tridiagonal sys
li=np.zeros(n-1)
di=np.zeros(n)
ui=np.zeros(n-1)
gamma=np.zeros(n)

m1 = (1/h**2)-(p/(2*h))
m2 = (-2/h**2)+q
m3 =  (1/h**2)+(p/(2*h))
for i in range(1,n-1):
  gamma[i]=np.e**(x[i+1])
  gamma[0]=np.e**(x[1])-m1*yL
  gamma[n-1]=np.e**(x[n])-2*m3*h*yR
# print(x[i+1])
print(gamma)
# gamma = np.exp(x)

# for i in range(0, n-1):
li[:] = m1
di[:] = m2
ui[:] = m3

# Sparse matrix A setup
A = diags([li, di, ui], [-1,0,1], shape=(n, n)).toarray()
A[n-1,n-2] =m1+m3

print(A)

# Solve the linear system
rslt = spsolve(A,gamma)
print('rslt=',rslt)


# Exact solution
u_exact= np.zeros(n)
for i in range(n):
    u_exact[i] = np.e**(x[i])+(1-np.e**(1))*x[i]

print('u_exact=',u_exact)

# Calculate and print the errors
errors = np.abs(rslt - u_exact)
for i in range(0,n):
    print(f"x = {x[i+1]:.1f}, Numerical = {rslt[i]:.4f}, Exact = {u_exact[i]:.4f}, Error = {errors[i]:.4e}")
# Plot the numerical and exact solutions
plt.rcParams["font.family"] = "times new roman"
plt.plot(x[:-1], rslt, label='Numerical Solution', marker='o')
plt.plot(x[:-1], u_exact, label='Exact Solution', marker='x')
plt.xlabel('x')
plt.ylabel('Solution')
plt.legend()
plt.title('Numerical vs Exact Solution for BVP')
plt.grid(True)
# plt.savefig('Q11',dpi=500)
plt.show()


#_________________________________________________________________________________________________________

# linear ODE constant Coeficient 
# both neumann BVP
# y′′+2y′+y=ex with y′(0)=1andy′(1)=1
# given parameters

import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
import math

a=0
b=1
yR = 1
yL = 1
n=5
# n = 4d
h=(b-a)/n
print(h)
p=2
q=1
# boundary value conditions
x=np.linspace(a,b,n+1)
print(x)
# tridiagonal sys
li=np.zeros(n)
di=np.zeros(n+1)
ui=np.zeros(n)
gamma=np.zeros(n+1)
print(gamma)
m1 = (1/h**2)-(p/(2*h))
m2 = (-2/h**2)+q
m3 =  (1/h**2)+(p/(2*h))
for i in range(1,n):
   gamma[i]=np.e**(x[i])
   gamma[0]=np.e**(x[0])+2*m1*h*yL
   gamma[n]=np.e**(x[n])-2*m3*h*yR

# for i in range(1,n):
   # gamma[i]=np.sin(x[i])
# gamma[0]=np.sin(x[0])+2*m1*h*yL
# gamma[n]=np.sin(x[n])-2*m3*h*yR
  
# print(x[i+1])
print(gamma)
# gamma = np.exp(x)

# for i in range(0, n-1):
li[:] = m1
di[:] = m2
ui[:] = m3

# Sparse matrix A setup
A = diags([li, di, ui], [-1,0,1], shape=(n+1, n+1)).toarray()
A[n,n-1] = m1+m3
A[0,1] = m1+m3    

print(A)

# Solve the linear system
rslt = spsolve(A,gamma)
print('rslt=',rslt)

# ___________________________________________________________________________________
# linear ODE constant Coeficient 
# right dirichlet and left neumann BVP
# y'' = e^x with y(0)=1 and y'(1)=1 Exact Solution e^x(1-e)x
# given parameters

a=0
b=1
yL = 1
yR =1 
n= 5
h=(b-a)/n
print(h)
p=0
q=0
# boundary value conditions
x=np.linspace(a,b,n+1)
print(x)
# tridiagonal sys
li=np.zeros(n-1)
di=np.zeros(n)
ui=np.zeros(n-1)
gamma=np.zeros(n)

m1 = (1/h**2)-(p/(2*h))
m2 = (-2/h**2)+q
m3 =  (1/h**2)+(p/(2*h))
for i in range(1,n-1):
  gamma[i]=np.e**(x[i])
  gamma[0]=np.e**(x[0])+2*m1*yL*h
  gamma[n-1]=np.e**(x[n-1])-m3*yR
# print(x[i+1])
print(gamma)
# gamma = np.exp(x)

# for i in range(0, n-1):
li[:] =m1
di[:] = m2
ui[:] = m3

# Sparse matrix A setup
A = diags([li, di, ui], [-1,0,1], shape=(n, n)).toarray()
A[0,1] =m1+m3

print(A)

# Solve the linear system
rslt = spsolve(A,gamma)
print('rslt=',rslt)


# Exact solution
u_exact= np.zeros(n)
for i in range(n):
    u_exact[i] = np.e**(x[i])+(1-np.e**(1))*x[i]

print('u_exact=',u_exact)

# Calculate and print the errors
errors = np.abs(rslt - u_exact)
for i in range(0,n):
    print(f"x = {x[i+1]:.1f}, Numerical = {rslt[i]:.4f}, Exact = {u_exact[i]:.4f}, Error = {errors[i]:.4e}")
# Plot the numerical and exact solutions
plt.rcParams["font.family"] = "times new roman"
plt.plot(x[:-1], rslt, label='Numerical Solution', marker='o')
plt.plot(x[:-1], u_exact, label='Exact Solution', marker='x')
plt.xlabel('x')
plt.ylabel('Solution')
plt.legend()
plt.title('Numerical vs Exact Solution for BVP')
plt.grid(True)
# plt.savefig('Q11',dpi=500)
plt.show()