import numpy as np
from math import erf

e=10**-3
k=np.sqrt(2*e)

def fun(x, u):
    return [u[1],-x*u[1]/float(e)-np.pi**2*np.cos(np.pi*x)-np.pi*x*np.sin(np.pi*x)/float(e)]

def bc(ua, ub):
    return [ua[0]+2, ub[0]]
def exact(x):
    return np.cos(np.pi*x)+erf(x/float(k))/float(erf(1./k))

x = np.linspace(-1, 1, 5)
u_a = np.zeros((2, x.size))

from scipy.integrate import solve_bvp
res_a = solve_bvp(fun, bc, x, u_a)
 
x_plot = np.linspace(-1, 1, 50)
u_plot_a = res_a.sol(x_plot)[0]

x = np.linspace(-1, 1, 50)
u_exact=np.zeros(len(x))
for i in range(len(x)):
    u_exact[i]=exact(x[i])

import matplotlib.pyplot as plt
plt.plot(x_plot, u_plot_a, 'bo-', label='u_a')
plt.plot(x, u_exact, 'g-', label='u_exact')
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.grid('on')
plt.savefig('Ex5.png')
plt.show()





