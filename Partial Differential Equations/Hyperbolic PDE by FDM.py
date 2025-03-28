import numpy as np 
import matplotlib.pyplot as plt
import time 
from scipy.sparse import diags
import scipy
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Wave Equation Explcit Method 
u_tt =c u_xx , c=1
u(x,0)=sin(πx) and ∂u/∂t(x,0)=0
u(0,t)=u(1,t)= 0 for all t.

h = 0.25
k = 0.1
c = 1
a = 0
b = 1

M = 4 # Number of spatial steps
N = 10 #time levels
x = np.linspace(a,b,M+1)
#lambda 
L = (c*k)/h

#inital condition
ui = np.sin(np.pi*x)
gx = 0*x
#boundary condition
ui[0] = 0
ui[M] = 0
# print(ui)
# explcit scheme
u = np.zeros([N+1,M+1])
# # print(u)
for j in range(0,N+1):
    for i in range(0,M+1):
        u[0,i] = ui[i]
        u[j,0] = 0
        u[j,M] = 0       

for i in range(1,M):          
    u[1,i] =( 2*u[0,i]+2*k*gx[i]+(L**2)*(u[0,i-1]-2*u[0,i]+u[0,i+1]))/2
    
for j in range(1,N):
    for i in range(1,M):
        u[j+1,i] = 2*u[j,i]-u[j-1,i]+L**2*(u[j,i-1]-2*u[j,i]+u[j,i+1])
        
print(u)

# # Plot the results
# plt.figure(figsize=(8, 6))
# for i in range(N+1):
#     plt.plot(x, u[i, :], label=f"t = {i*k}", marker = 'o')

# plt.xlabel('x')
# plt.ylabel('u(x,t)')
# plt.title('Explicit Scheme Solution to the Heat Equation')
# plt.legend(loc='upper right')
# plt.grid(True)
# plt.show()

# # Enhanced plot with Times New Roman font and improved quality
# plt.rcParams.update({
#     "font.family": "Times New Roman",  # Set font to Times New Roman
#     "font.size": 12,                   # Set font size
#     "figure.dpi": 300                  # Higher resolution for better quality
# })

plt.figure(figsize=(8, 6))
for i in range(N + 1):
    plt.plot(x, u[i, :], label=f"t = {i}", marker='o')  # Use integer time labels

plt.xlabel('x')
plt.ylabel('u(x,t)')
plt.title('Explicit Scheme Solution to the Wave Equation')
plt.legend(loc='upper right')
plt.grid(True)
plt.show()

