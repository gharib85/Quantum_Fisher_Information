
#Assumptionless SineState
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
                s_coeff = (s_coeff+(-1)**(n-k+s)*(np.sqrt(float(scipy.special.factorial(n,exact=True)))
                                                  /scipy.special.factorial(n-k+s,exact=True))*
                (np.sqrt(float(scipy.special.factorial(N-n,exact=True)))/(scipy.special.factorial(N-n-s,exact=True)))*
                (np.sqrt(float(scipy.special.factorial(k,exact=True)))/(scipy.special.factorial(k-s,exact=True)))*
                (np.sqrt(float(scipy.special.factorial(N-k,exact=True)))/(scipy.special.factorial(s,exact=True))))
            
        state[n]=(state[n]+s_coeff*(np.sin(np.pi/4.)**N)*
                      np.sin((k+1)*np.pi/(N+2))*np.e**(I*np.pi*(k-n)/2))
    
    state[n]=state[n]/np.sqrt(float(N*0.5+1))

#print("State = %s" % np.round(state,4))
SCS = np.real((state*np.conj(state)))

state_expanded = []
for t in xrange(0,N+1): 
    state_expanded = np.append(state_expanded, np.ones(int(scipy.special.binom(N, t)))*
                                                            (state[t]/np.sqrt(float(scipy.special.binom(N, t)))))
    
    #state_expanded = np.round(state_expanded,4)
    
state_expanded_mod = np.real((state_expanded*np.conj(state_expanded)))
for h in xrange(1,2):

    r=0.05
    t=1.-r

#plt.plot(list(range(0,N+1)), np.round(SCS,4), 'ro')
#plt.xlabel(r'$ \nu $')
#plt.ylabel('Prob')
#plt.title('N = %s and r = %s (Input)' %(N,r))
#plt.show()

    #Strings = range(N+1)
#Output_Prob = np.zeros((N+1,N+1))
    Hamming = np.zeros(N+1)

    #for i in xrange(0,N+1):
        #Strings[i]= [int(j) for j in list('0'*(N-i)+'1'*i)]

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

    for i in xrange(0,2**N):
        for j in xrange(0,2**N):
            Trans = len(compare_listcomp(comb[i], comb[j]))
            Hamming[int(sum(comb[j]))] = (Hamming[int(sum(comb[j]))] + 
            (t**Trans)*(r**(N-Trans))*state_expanded_mod[i])
    
    Hamming = np.divide(Hamming, np.sum(Hamming))
    
#print(Hamming)
    plt.plot(list(range(0,N+1)), np.round(SCS,4), 'ro', label='Sine State')
    plt.plot(list(range(0,N+1)), Hamming, 'bo', label='Hamming')
    plt.legend(loc='upper left', frameon=False)
    plt.xlabel('Units')
    plt.ylabel('Prob')
    plt.title('N = %s and r = %s' %(N,r))
    plt.show()
