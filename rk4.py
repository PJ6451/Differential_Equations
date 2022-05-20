import numpy as np
import matplotlib.pyplot as plt

def rk4(func, x, dt): 
    k1 = dt*func(x)
    k2 = dt*func(x + k1*0.5)
    k3 = dt*func(x + k2*0.5)
    k4 = dt*func(x + k3)
    return x + (1/6.)*(k1 + 2.*k2 + 2.*k3 + k4)

def numsolv(func, x0, dt, Niter):
    soln = np.zeros((Niter,1),dtype=float)
    soln[:,0]
    for jj in range(Niter):
        x0 = rk4(func,x0,dt)
        soln[jj] = x0
    return soln
