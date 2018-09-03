
# coding: utf-8

# In[20]:


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
#Scaling = np.zeros(MAX-1)
for d in xrange(1,theta_max+1):
    for N in xrange(1,MAX+1):
        state = np.array([complex(item) for item in np.zeros(N+1)])

        for n in xrange(N+1):
    
            for k in xrange(N+1):
        
                if (n-k)<0:
                    s_min=abs(n-k)
                else:
                    s_min=0     
                if k > (N-n):
                    s_max=N-n
                elif k <= (N-n):
                    s_max=k
            
                s_coeff=0
                for s in xrange(s_min,(s_max+1)):
                    s_coeff = (s_coeff+(-1)**(n-k+s)*(np.sqrt(float(scipy.special.factorial(n,exact=True)))
                                                  /scipy.special.factorial(n-k+s,exact=True))*
                    (np.sqrt(float(scipy.special.factorial(N-n,exact=True)))/(scipy.special.factorial(N-n-s,exact=True)))*
                    (np.sqrt(float(scipy.special.factorial(k,exact=True)))/(scipy.special.factorial(k-s,exact=True)))*
                    (np.sqrt(float(scipy.special.factorial(N-k,exact=True)))/(scipy.special.factorial(s,exact=True))))
            
                state[n]=state[n]+s_coeff*(np.sin(np.pi/4.)**N)*np.sin((k+1)*np.pi/(N+2))*np.e**(I*np.pi*(k-n)/2)
    
            state[n]=state[n]/np.sqrt(float(N*0.5+1))

        #print("State = %s" % np.round(state,4))
        #print("Total prob. = %s" % round(np.sum(np.real((state*np.conj(state)))),4))

        state_expanded = []
        for t in xrange(0,N+1): 
            state_expanded = np.append(state_expanded, np.ones(int(scipy.special.binom(N, t)))*
                                                               (state[t]/np.sqrt(float(scipy.special.binom(N, t)))))
    
        state_expanded = np.round(state_expanded,4)
        #print ("State in full = %s" % np.round(state_expanded,4))

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
        OS=Operator*(np.matrix(state_expanded).T)
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
