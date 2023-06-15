import random
from numpy import *
from numpy.fft import*
Z= [1,0,1,0,1,0,1,0,1,0,1,1,1,0]
Y=ifft(Z)
print(Y[0])
#print(Y)
for i in range(len(Y)):
    Z[i] = abs(Y[i]**2)
#print(Z)
Z = Z/sum(Z)
print(sum(Z))
#print(Z)
p1 = 29
p2 = 31
N  = p1 * p2 # N não precisa ser semi-primo
x  = 1024
r  = 0
q  = 2**20
#result=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
#K=EstimaOrdem(r, result)


