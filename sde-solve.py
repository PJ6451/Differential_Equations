import numpy as np
import sdeint
import matplotlib.pyplot as plt

#Parameters
a = 1.0
b = 1.0
s = 1.0

#Conditions and timespan
tspan = np.linspace(0,10,1000)
x0 = 0.1

#defining functions
def f(x,t):
    return a*(b-x)

def g(x,t):
    return s

#numerical solver
result = (sdeint.itoint(f,g,x0,tspan))

#plot
plt.plot(tspan,result)