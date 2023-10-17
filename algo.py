import random
from numpy import *
from numpy.fft import *
import time
import math
import sys
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
    return L
def Frac(L) :
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
            #print(F[1][2:m])
            F[1][0]=lcm.reduce(array(F[1][2:m]))
            print('Fzin',F[1][0])
    if F[1][1:m]==[]:
        print('faz parte')
        #print("Fração inexistente:",L)

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
            print('result ',result[i]+j)
            print(FracCont(result[i]+j,q))
            t=Frac(FracCont(result[i]+j,q))
            l[j+1]=t[1][0]
            print('tzin: ',t[1][0])
            print(sys.getsizeof(t[1][0]))
        R.append(l)
    return R

p1 = 23
p2 = 29
N  = p1 * p2 # N não precisa ser semi-primo
x  = 1028150
r  = 308
q  = 2**20
result= [1, 74898, 779624, 953252, 1031555, 289380, 1, 105539, 422155, 255336, 1021340, 674085, 1, 1, 735365]

R= EstimaOrdem(r,result)
print(R)