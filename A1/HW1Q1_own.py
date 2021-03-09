# want to find min cum energy
# given structure s and sequences w1 w2 w3...
#

# want s st. min of (max energy of s for 3 ws)
import numpy as np
import math

def score(i,j, seq, struct):
    #check if should bp closed parenthesis
    if(struct[i]!= '(' or struct[j]!=')'):
        return 0
    if struct[i:j+1].count('(') != struct[i:j+1].count(')'):
        return 0
    stack = []
    for val in struct[i:j+1]:
        if val == '(':
            stack.append('(')
        elif val == ')':
            try:
                stack.pop()
            except:
                return 0
    if stack: #not empty stack
        return 0

    # if should bp
    if (seq[i], seq[j]) in [('A','U'), ('U', 'A')]:
        return -2
    elif (seq[i], seq[j]) in [('G', 'C'), ('C', 'G')]:
        return -4
    elif (seq[i], seq[j]) in [('G', 'U'), ('U', 'G')]:
        return -1
    else: # if bp and not match, is incompatible
        print('INFINITY"')
        return np.inf # since want to minimize
        #NOT for s1 or s2 (for s3)

def getLs(struct, omegas):
    # omegas = [omegas[0]]
    # struct = '(.)'
    # omegas = ['ACU']
    print(omegas)
    print(struct)

    #set M
    theta = 1
    N = len(struct)
    M = np.empty((N,N))
    M[:] = np.NAN
    for k in range(0, theta):
        for i in range(N-k):
            j = i + k
            M[i][j] = 0 # set diagonal to 0

    #get L(s) energy
    for L in range(theta, N):
        for i in range(0, N - L):
            j = i+L
            if L <= theta:
                M[i][j] = 0
            #no bp at j
            maxval = M[i][j-1] # if j unpaired
            for k in range(i, j-theta):
                print(i,k, j, 'K')
                maxscore = -np.inf
                for omega in omegas:
                    # max b/c want least negative
                    maxscore = max(maxscore, score(k,j,omega, struct))
                    print(maxscore)
                print(maxscore, 'MAX')
                maxval = max(maxval, np.abs(maxscore))
                # print(newE, 'E')
                # print(newE)
                print(maxval, 'max', maxval, maxscore, M[i][k-1] + M[k+1][j-1] + maxscore)
                M[i][j] = maxval
            # print(M[i][j])
            #     print('new', M[i][j])
            #     print(i, j, M[i][j])
            # print(M[i][j])
        # if L > 2:
        #     break
    print(M)


def main():
    # structs = ['((..))((..))', '((((....))))', '((..))((..))']
    structs = ['((..))((..))']
    omegas = ['GGAAUUGGUUCC', 'GGAAUU-GUUCC', 'AAAAUUGGUUUU']
    for struct in structs:
        getLs(struct, omegas)
        break

if __name__ == '__main__':
    main()
