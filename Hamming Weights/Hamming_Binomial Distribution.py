
# coding: utf-8

# In[2]:


from sympy import Matrix, I
from sympy import *
import numpy as np
from sympy.physics.quantum import TensorProduct
import matplotlib.pyplot as plt
from matplotlib import pyplot
#import math
import scipy.special
from scipy import stats
#import itertools
#from itertools import combinations

#N = input("N = ")
N=15
state_before = range(N+1)
for i in xrange(0,N+1):
    state_before[i] = scipy.special.binom(N, i)/2**N

for h in xrange(0,20):
    r=h*0.05
    t=1.-r

    Strings = range(N+1)
#Output_Prob = np.zeros((N+1,N+1))
    Hamming = np.zeros(N+1)

    for i in xrange(0,N+1):
        Strings[i]= [int(j) for j in list('0'*(N-i)+'1'*i)]

#print(Strings)

    def generate_binary(n, l):

      if n == 0:
        return l
      else:
        if len(l) == 0:
          return generate_binary(n-1, ["0", "1"])
        else:
          return generate_binary(n-1, [i + "0" for i in l] + [i + "1" for i in l])

    comb = [list(j) for j in generate_binary(N , [])]

    for i in xrange(0,2**N):
        comb[i]= [int(j) for j in comb[i]]

#print(comb)

    def compare_listcomp(x, y):
        return [i for i, j in zip(x, y) if i == j]

    for i in xrange(0,N+1):
        for j in xrange(0,2**N):
            Trans = len(compare_listcomp(Strings[i], comb[j]))
            Hamming[int(sum(comb[j]))] = (Hamming[int(sum(comb[j]))] + 
            (t**Trans)*(r**(N-Trans))*scipy.special.binom(N, i)/2**N)

#print(Hamming)
    plt.plot(list(range(0,N+1)), state_before, 'ro', label='Binomial')
    plt.plot(list(range(0,N+1)), Hamming, 'bo', label='Hamming')
    plt.legend(loc='upper left', frameon=False)
    plt.xlabel('Units')
    plt.ylabel('Prob')
    plt.title('N = %s and r = %s' %(N,r))
    plt.show()
    
    print ("Cumulants (Hamming): %s, %s, %s, %s" % (stats.kstat(Hamming, 1), stats.kstat(Hamming, 2), 
       stats.kstat(Hamming, 3), stats.kstat(Hamming, 4)))
