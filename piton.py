import random
from numpy import *
from numpy.fft import *
import time

def FracCont(x,q):
    x = abs(x)
    x_inic = x
    if x == 1 :
        L=[1,0]
        x = 0
    else:
        L=[]
    x = x/q
    i=0
    max = log(x_inic)/log(1.6)  # baseado nos termos da série de Fibonacci, para estimar a quantidade máxima de termos da fração continuada.
    print('max: ',max)
    while x>0  and i< max:
        c = int(x)
        L.append(c)
        x -= c
        if x >=0 :
           x = 1/x
           if x>x_inic or x>N:
               x=0
        else :
           x=0
        i += 1
    print('iiii',i)
    return L

# Recebe uma fração contínua e monta as frações aproximadas
'''def Frac(L) :
    m = len(L)
    if m==0:    # leu 1
        m=1
        L=[1]
    F=[[0]*m,[0]*m]
    num=1
    den=1
    for i in range(m-1):
        num=1
        den=L[m-i-1]
        for j in range(m-i-2,0,-1):
            num, den = den, L[j]*den+num
        F[0][i+1]=num
        F[1][i+1]=den
    F[0][0]=m-1
    if F[1][2:m]!=[]:
        if m>1:
            F[1][0]=lcm.reduce(array(F[1][2:m]))
    if F[1][1:m]==[]:
            print("Fração inexistente:",L)

    print('\n\ntotal: ',F[1][0], '\n\n')
    return F

def EstimaOrdem(r,result):
    mult = 0
    n    = len(result)
    R    = []

    print('n:',n)
    print('Tenta estimar a ordem r=ord(x,N) ou múltiplo ou divisor dela para extrair os fatores de N')
    for i in range(n):
        l=[1,1,1]
        for j in range(-1,2):
            print('.',end='')
            t=Frac(FracCont(result[i]+j,q))
            l[j+1]=t[1][0]
        R.append(l)
    return R'''
# #------------------------- Início do programa ------------------------
p1 = 47
p2 = 43
N  = p1 * p2 # N não precisa ser semi-primo
x  = 2
r  = 322
q  = 2**20
n  = 15 # quantidade de valores medidos
result=[26051, 179105, 696880, 898779, 113976, 709905, 1045319, 719675, 436363, 169335, 280054, 856446, 179105, 807599, 595929]
l =FracCont(x,q)
print(l)
#R=EstimaOrdem(r,result)
#print('\nvalor de R: \n', R)
