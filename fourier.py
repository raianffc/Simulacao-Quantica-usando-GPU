import random
from numpy import *
from numpy.fft import *

Z = [0, 1]
Y=ifft(Z)
print(Y)
print(Z)