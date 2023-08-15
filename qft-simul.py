import random
from numpy import *
from numpy.fft import *
import time

def buscabin(Soma,P,m):
    n = len(Soma)
    if n==0 :
        return 0
    elif n==1 :
        return P[0]
    elif n==2 :
        if m <= Soma[0]:
            return P[0]
        else :
            return P[1]
    else :
        meio = n//2
        if m == Soma[meio] :
            return P[meio]
        elif m < Soma[meio] :
            return buscabin(Soma[:meio],P[:meio],m)
        else :
            return buscabin(Soma[meio:],P[meio:],m)

# usar apenas para valores positivos
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
    if F[1][2:m]!=[]:
        if m>1:
            F[1][0]=lcm.reduce(array(F[1][2:m]))
    if F[1][1:m]==[]:
            print("Fração inexistente:",L)

    return F


#
# def mostraQFT(T) :
# # mostra o gráfico, caso seja desejável
#     clf(Z)
#     plot(Z)
#

def Prepara(N,x,r=0,q=1):
# colocar o valor de r caso seja conhecido (evita o cálculo abaixo)
    tamN = int(log2(N))
    q1 = 2**(2*tamN)    # este é o valor ideal segundo Shor. Não é usado no programa. Serve apenas de referência
    print("valor ideal para q:",q1)
    if(q < N) :
        q = 1 << (tamN+4)
    if r==0 :
        s=x
        i=1
        while s > 1 :
            s = s*x % N
            i += 1
        r = i
        print('Ordem r não informada. Ordem r calculada:',r)
    else:
        print('Ordem r informada:',r)
    print('Criando Z...')
    Z = [0]*q
    j = 1
    cont =0
    while  j <= q  :
        Z[j] = 1
        j += r
        cont += 1  # em princípio, não é usado. Mas poderia ser usado para normalizar o vetor de estado
#    print(cont)
    print('Calculando FFT...')
    Y=ifft(Z)
    print(Y)
    print('Calculando probabilidades...')
    for i in range(len(Y)):
        Z[i] = abs(Y[i]**2)
    Z = Z/sum(Z)    # normaliza a probabilidade
    print(sum(Z))
    print(Z)
#    mostraQFT(Z)

    print('Criando Soma com probabilidade acumulada')
    P    = [0]*(2*r+1)
    Soma = [0]*(2*r+1)
    k    = q/r
    print('Calculando Somas...')
    total = 0
    for i in range(r):
        pos= int(i*k)
        P[2*i] = pos
        total += Z[pos]
        Soma[2*i]= total
        P[2*i+1]= pos+1
        total = total + Z[pos+1]
        Soma[2*i+1] = total
    P[2*r-1]=-int(random.random()*q)
    Soma[2*r-1]=1
    print('Probabilidade Acumulada (dos picos):',total)

    return [Soma,P,r]

def Simula(Soma,P,n):
    print('Simula medição:')
    result = [0]*n
    for i in range(n):
        m=random.random()
        result[i] = buscabin(Soma,P,m)
        if result[i]==0:
            result[i] = 1

    return result

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

# 0 - não é múltiplo e nem divisor da ordem que distinga fator
# 1 - múltiplo da ordem
# 2 - distingue fator
# 4 - mmc múltiplo da ordem (acumulada)
# 8 - mmc distingue fator
def EstimaFator(N,x,R):
    print('Procura um múltiplo da ordem ou um divisor que distiga um fator não trivial.')
    Sucesso=[]
    potTotal=1
    for i in range(len(R)):
        sucesso=[0,0,0]
        for j in range(3):
            pot=pow(x,int(R[i][j]),N)
            if pot==1:
                sucesso[j]=1
                print('múltiplo de r')
            else:
                d=gcd(pot-1,N)
                if d > 1:
                    sucesso[j]=2
                    print('fator:',d)
            potTotal *= pot
            if potTotal % N ==1:
                sucesso[j]+=4
                print('fatores do múltiplo da ordem',pot,potTotal//pot)
            else:
                d=gcd(potTotal-1,N)
                if d > 1:
                    sucesso[j]+=8
            potTotal %= N
            # print(potTotal)
        Sucesso.append(sucesso)
    return Sucesso


def Fatores(N,x,R,S):
    fat=set()
    for i in range(len(S)):
        for j in range(3):
            if S[i][j]==1: # o valor é um múltiplo da ordem
                print('testa um divisor da ordem')
                if R[i][j] % 2 == 0:
                    print('Teste de Shor...')
                    f = gcd(pow(x, int(R[i][j] // 2), N) - 1, N)
                    fat.add(f)
                    fat.add(N // f)
                if R[i][j] % 3 == 0:
                    print('Teste com 3 ...')
                    f = gcd(pow(x, int(R[i][j] // 3), N) - 1, N)
                    fat.add(f)
                    fat.add(N // f)
            elif S[i][j] == 2:
                f = gcd(pow(x, int(R[i][j]), N) - 1, N)
                fat.add(f)
                fat.add(N // f)
    return fat
# #------------------------- Início do programa ------------------------
inicio = time.time()
p1 = 47
p2 = 43
N  = p1 * p2 # N não precisa ser semi-primo
x  = 2
r  = 0
q  = 2**20
Soma,P,r=Prepara(N,x,r,q)
#print("sum",Soma)
n  = 15 # quantidade de valores medidos
result=Simula(Soma,P,n)
print('q=',q)
print('Resultados medidos na rotina quântica do QOFA:')
print(result)
R=EstimaOrdem(r,result)

S=EstimaFator(N,x,R)
fat=Fatores(N,x,R,S)
print(fat)
fim = time.time()
print('Tempo de execução: ', fim-inicio)