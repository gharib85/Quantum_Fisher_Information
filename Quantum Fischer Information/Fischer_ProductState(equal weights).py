
# coding: utf-8

# In[21]:


from sympy import Matrix, I
from sympy import *
import numpy as np
from sympy.physics.quantum import TensorProduct
#import matplotlib.pyplot as plt
#import math
from scipy.special import factorial
import scipy.special

MAX=6
theta_max = 12 
FischerNum = np.zeros(MAX)
#Scaling = np.zeros(MAX)
for d in xrange(1,theta_max+1):
    for N in xrange(1,MAX+1):
        
        state = np.divide(np.ones(2**N),2**(N*0.5))
        #print("State = %s" % np.round(state,4))
        #print("Total prob. = %s" % round(np.sum(np.real((state*np.conj(state)))),4))

        PHI = np.pi/3.0
        THETA = np.pi*d/26.0
        theta = Symbol('theta', real = True)
        Operator0 = Matrix([[cos(PHI) - I*cos(2*theta)*sin(PHI), -sin(2*theta)*sin(PHI)], [sin(2*theta)*sin(PHI), 
                                                                            cos(PHI) + I*cos(2*theta)*sin(PHI)]])

        s1="TensorProduct("
        for i in xrange(N):
            s1=s1+"Operator0,"
        s1=s1+")"
        Operator = eval(s1)
        OS=Operator*(np.matrix(state).T)
        OSC=OS.C
        Prob = np.multiply(OS,OSC)
        #print ("Total prob. after applying the operator = %s" % 
               #complex((simplify(np.sum(Prob).subs({theta: THETA})))).real)

        Fischer = 0
        for i in range(0, 2**N):
            Fischer = Fischer + (OS[i]*OSC[i])*(((log((OS[i]*OSC[i]))).diff(theta))**2)
        FischerNum[N-1] = complex(Fischer.subs({theta: THETA})).real
        #print ("Fischer Information = %s (phi = %s, theta = %s)" % (complex(FischerNum).real, PHI, THETA))
        #Scaling[N-1] = (FischerNum[N-1])**-0.5
        #print ("Del theta = %s" % Scaling)

        #print ("FischerNum(N = %s) = %s" % (N,FischerNum[N-1]))
    r = (np.cos(THETA/2))**2
    print ("FischerNum (theta = %s pi, r = %s) = %s" % (THETA/np.pi, r, FischerNum))
    #print ("Scaling = %s" % Scaling)
