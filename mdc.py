import random
from numpy import *
from numpy.fft import *

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

# Recebe uma fração contínua e monta as frações aproximadas
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
    print('print F[0][0]',F[0][0])
    if F[1][2:m]!=[]:
        if m>1:
            F[1][0]=lcm.reduce(array(F[1][2:m]))
    if F[1][1:m]==[]:
            print("Fração inexistente:",L)

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
    return R
# #------------------------- Início do programa ------------------------
p1 = 37
p2 = 41
N  = p1 * p2 # N não precisa ser semi-primo
x  = 2
r  = 0
q  = 2**20

result = [343700, 774782, 570891, 367001, 93206, 1042750, 425256, 1, 757304, 1042750, 343700, 1042750, 297096, 75730, 413604]
R = EstimaOrdem(r, result)
print(R)