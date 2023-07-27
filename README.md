# Simulação de Computadores Quânticos usando GPU
Estou traduzindo o código em Python para linguagem C. Depois irei paralelizar tudo em CUDA.
Este é um código que implementa o algoritmo de Shor para fatoração de números inteiros. Consiste em várias funções que são usadas para executar as diferentes etapas do algoritmo.

buscabin(Soma,P,m): Esta função realiza uma busca binária em uma lista de probabilidades para encontrar o índice do primeiro elemento maior ou igual a um determinado número m. É usado para simular uma medição de um estado quântico no algoritmo.

FracCont(x,q): Esta função calcula a expansão contínua da fração de um número racional x em relação a um determinado denominador q. Ele retorna uma lista dos quocientes parciais da fração contínua.

Frac(L): Esta função pega uma lista de quocientes parciais de uma fração contínua e retorna uma lista dos convergentes correspondentes (frações que se aproximam do número racional original).

Prepara(N,x,r=0,q=1): Esta função prepara o estado quântico usado no algoritmo de Shor. Ele toma como entrada o inteiro a ser fatorado N, um inteiro aleatório x que é primo de N, a ordem r de x módulo N (se conhecido) e um parâmetro opcional q que determina o número de elementos no estado quântico. Ele retorna uma lista de duas listas: Soma (uma lista das probabilidades acumuladas do estado quântico) e P (uma lista dos índices dos elementos do estado quântico).

Simula(Soma,P,n): Esta função simula a medição do estado quântico n vezes e retorna uma lista dos resultados.

EstimaOrdem(r,result): Esta função pega a ordem r de x módulo N e uma lista de resultados de medição e os usa para estimar os fatores de N.

O algoritmo principal não está incluído no código e precisaria ser implementado separadamente. Isso envolveria chamar essas funções na ordem apropriada para executar a parte quântica do algoritmo, seguida pelo pós-processamento clássico dos resultados da medição para obter os fatores de N.

Sites sobre a biblioteca Complex.h 

https://learn.microsoft.com/pt-br/cpp/standard-library/complex?view=msvc-170

https://acervolima.com/arquivo-de-cabecalho-complex-h-em-c-com-exemplos/

https://www.geeksforgeeks.org/complex-h-header-file-in-c-with-examples/

https://www.inf.ufpr.br/roberto/ci067/18_casting.html

https://www.codeproject.com/Articles/9388/How-to-implement-the-FFT-algorithm

https://www.katjaas.nl/FFTimplement/FFTimplement.html

Agora para que o código seja mais otimizado usaremos o FFTW uma biblioteca feito por pesquisadores do MIT.

https://www.fftw.org/

Biblioteca ja instalada, vendo como usar IFFT.
Tive problemas com algoritimo em sim. por termos uma lista, mesmo que unidimensional, muito grande e demorava muito para executar. Perdendo o próposito do projeto.

Normalizar o vetor esta sendo um problema. Não consigo tirar elevar o quadrado para tirar o conjugado.
