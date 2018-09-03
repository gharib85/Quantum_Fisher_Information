
# coding: utf-8

# In[1]:


from sympy import Matrix, I
from sympy import *
import numpy as np
from sympy.physics.quantum import TensorProduct
import matplotlib.pyplot as plt
from matplotlib import pyplot
import math
import random
import scipy.special
from scipy import stats

N = input("N = ")
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
            s_coeff = (s_coeff+(-1)**(n-k+s)*(math.sqrt(math.factorial(n))/math.factorial(n-k+s))*
            (math.sqrt(math.factorial(N-n))/(math.factorial(N-n-s)))*
            (math.sqrt(math.factorial(k))/(math.factorial(k-s)))*
            (math.sqrt(math.factorial(N-k))/(math.factorial(s))))
            
        state[n]=state[n]+s_coeff*(math.sin(math.pi/4)**N)*math.sin((k+1)*math.pi/(N+2))*math.e**(I*math.pi*(k-n)/2)
    
    state[n]=state[n]/math.sqrt(N*0.5+1)

#print("State = %s" % np.round(state,4))
#print("Total prob. = %s" % round(np.sum(np.real((state*np.conj(state)))),4))

state_expanded = []
for t in xrange(0,N+1): 
    state_expanded = np.append(state_expanded, np.ones(int(scipy.special.binom(N, t)))*
                                                               (state[t]/math.sqrt(scipy.special.binom(N, t))))
    
state_expanded = np.round(state_expanded,4)

print("State_expanded = %s" % state_expanded)
