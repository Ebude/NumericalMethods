import numpy as np
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve, norm

def Bvp(fcn, a, b, ua, ub, Nx):
    """ 
    Solution of 1D boundary value problem 
        -u" = fcn(x) for -1 < x < 1
    with boundary conditions
        u(-1) = 2,  u(1) = 0.
    on a uniform mesh
    """

    L = b-a                      # length of the domain
    dx = float(L)/float(Nx)      # length of the grid cell
    x = np.linspace(a, b, Nx+1)  # the grid
    u = np.zeros(Nx+1)           # the solution

    # The sparse matrix
    dia  = 2*np.ones(Nx-1) / dx**2
    low  = -np.ones(Nx-1) / dx**2
    upp  = -np.ones(Nx-1) / dx**2
    A = spdiags([low,dia,upp],[-1,0,1],Nx-1,Nx-1) 
    # print A.todense()

    # evaluate right hand side
    rhs = np.zeros(Nx-1)         # the right hand side
    for j in range(Nx-1):
        rhs[j] = fcn(x[j+1])
    rhs[0]  -= low[0]*ua
    rhs[-1] -= upp[0]*ub

    # solve the linear system
    u[0], u[-1]  = ua, ub
    u[1:-1] = spsolve(A, rhs)
    
    # return the grid and the numerical solution
    return u, x

if __name__ == '__main__':
    Nx = 20    # Nx = number of grid cells, 
                # Nx+1 = number of grid points
    a = -1.     # a = left end of the domain
    b = +1.     # b = right end of the domain
    ua =2     # boundary value left side
    ub = 0     # boundary value right side 

    def fcn(x):
        return 6*x

    def exact(x):
        return -x**3+1
   
    import time 
    t0 = time.clock()
    u,x = Bvp(fcn, a, b, ua, ub, Nx)
    t1 = time.clock()
    print 'cpu time',t1-t0 

    # compute the error norm
    from scipy.linalg import norm
    print 'approximation error', abs(norm(u-exact(x)))

    xx = np.linspace(a,b)
    sol = exact(xx)
 
    # generate a plot of the solution
    import matplotlib.pyplot as plt
    plt.plot(x,u,'bo-', linewidth=.5,label='Num')
    plt.plot(xx,sol,'r-',label='Exact')
    plt.legend()
    plt.title('numerical solution')
    plt.grid('on')
    plt.savefig('Ex4.png')
    plt.show()
    

