import numpy as np
import matplotlib.pyplot as plt

##Exercise 3b##
def fcn(u,t):
    return -3*u + 2*t

def exact(t): #solution of the ODE u'=f(u,t) using separation of variables
    return (2./3)*(t-(1./3))-(2./3)*(np.pi-(1./3))*np.exp(-3*(t-np.pi))
   
t0 = np.pi; te = 0 #initial and final time
Nt = 40
I = exact(0) #initial value of u

u = np.zeros(Nt+1)  
u1 = np.zeros(Nt+1)  
t = np.linspace(t0, te, Nt+1)
dt = t[0] - t[1] #stepsize

u[0] = I
u1[0]=I
for n in range(Nt):
    u[n+1] =u[n] + dt*fcn( u[n], t[n] ) #Euler's Method
    k1 = dt*fcn( u1[n], t[n] )          
    k2 = dt*fcn( u1[n]+k1, t[n+1] )
    u1[n+1] = u1[n] + 0.5 * ( k1 + k2 ) #Heun's method

tt = np.linspace(t0, te, 101)
ex = exact(tt)

plt.plot(t0-t, u,'o-', tt, ex,'r-',t0-t,u1,'g-')
plt.xlabel('t'); plt.ylabel('u')
plt.legend(['Euler','exact','Heun'], loc='lower left')
plt.title('Euler and Heun method, numerical solution')
plt.grid('on')
plt.savefig('Euler0.pdf',bbox_inches='tight')
plt.show()


##Exercise 4##

#a
aVec=np.array([3.14,15,9,26])
aVec

#b
bVec=np.arange(5,-5.2,-0.2)
bVec

##Exercise 5##
#a
aMat=2*np.ones((9,9))
aMat

#b
bMat=2*np.zeros((9,9))
bMat=np.diag([1,2,3,4,5,4,3,2,1])
bMat

#c
cMat=np.arange(100).reshape((10,10))
cMat=cMat.transpose()
cMat

#d
dMat=np.array([[13,-1,5],[-22,10,-87]])
dMat


