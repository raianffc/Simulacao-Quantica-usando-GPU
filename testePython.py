import random
from numpy import *
from numpy.fft import *

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
    #print(Y)
    print('Calculando probabilidades...')
    for i in range(len(Y)):
        Z[i] = abs(Y[i]**2)
        print(Z[i])
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

p1 = 31
p2 = 29
N  = p1 * p2 # N não precisa ser semi-primo
x  = 2
r  = 0
q  = 2**20

Soma,P,r=Prepara(N,x,r,q)
#print("sum",Soma)
n  = 15 # quantidade de valores medidos
'''
result=Simula(Soma,P,n)
print('q=',q)
print('Resultados medidos na rotina quântica do QOFA:')
print(result)

R=EstimaOrdem(r,result)

S=EstimaFator(N,x,R)
print('tamanho de S',len(S))
fat=Fatores(N,x,R,S)
print(fat)
'''