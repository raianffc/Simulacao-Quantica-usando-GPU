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


A biblioteca FFTW (Fastest Fourier Transform in the West) é uma biblioteca de transformada de Fourier altamente otimizada e amplamente usada. Ela oferece funções para calcular transformadas de Fourier diretas e inversas em várias dimensões e suporta diferentes tipos de dados. Aqui estão alguns dos comandos mais comuns usados na FFTW3:

Criação de um Plano:

fftw_plan fftw_plan_dft_1d(int n, fftw_complex *in, fftw_complex *out, int sign, unsigned flags);
fftw_plan fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags);
fftw_plan fftw_plan_dft_c2r_1d(int n, fftw_complex *in, double *out, unsigned flags);
Essas funções são usadas para criar um plano de transformada. n é o tamanho da transformada, in é o array de entrada, out é o array de saída (para transformadas complexas) ou NULL para transformadas reais, sign é 1 para FFT e -1 para IFFT, e flags são as opções de otimização.

Execução da Transformada:

void fftw_execute(const fftw_plan plan);
Essa função executa a transformada especificada pelo plano.

Liberando o Plano:

void fftw_destroy_plan(fftw_plan plan);
Essa função libera a memória alocada para um plano.

Alocação e Liberação de Memória:

fftw_complex* fftw_malloc(size_t size);
void fftw_free(void *ptr);
Use fftw_malloc para alocar memória para arrays de tipo fftw_complex (complexos) e fftw_free para liberar essa memória.

Configurações de Otimização:

#define FFTW_ESTIMATE 1
#define FFTW_MEASURE 2
#define FFTW_PATIENT 3
#define FFTW_EXHAUSTIVE 4
Essas macros são usadas para especificar os níveis de otimização ao criar um plano. FFTW_ESTIMATE é mais rápido, enquanto FFTW_EXHAUSTIVE é mais preciso.

para compilar no gcc no ubuntu 
-lfftw3 -lm -I/usr/local/lib

Agora terminar a parte de CUDA

Esses são apenas alguns dos comandos mais comuns usados na biblioteca FFTW3. A documentação oficial do FFTW3 fornece informações detalhadas sobre todas as funções disponíveis e opções de configuração. Certifique-se de verificar a documentação para obter mais informações: http://www.fftw.org/fftw3_doc/index.html

Para medir o tempo de execução do programa em C.
https://www.techiedelight.com/pt/find-execution-time-c-program/
Consertei o código em python, na função Fatores, na hora de fazer o mdc para saber se é ou nao fator de N, ele apenas printava e nao armazenava dentro da lista. Ultima coisa que poderia fazer para poder ficar mais rápido seria passar as listas para numpy para poder ficar mais rápido.
Nao é compativel com modelos GTX, apenas RTX 

Primeiro teste
  N = 29*31;
  "qbits" = 3;
  chute =2;
  q = 2^24;
Tempo de execução:
  Cuda: 3.16 segundos
  C: 0.8967 segundos
  Python: 6.0265 segundos

Segundo teste
  N = 29*31;
  "qbits" = 3;
  chute =2;
  q = 2^20;
Tempo de execução:
  Cuda: 0.159 segundos
  C: 0.0407 segundos
  Python: 0.3918 segundos

Terceiro teste
  N = 83*89;
  "qbits" = 3;
  chute =2;
  q = 2^24;
Tempo de execução:
  Cuda: 0.729 segundos Nao acertou, chute talvez
  C: 0.958328 segundos
  Python: 6.0848 segundos

Quarto teste
  N = 83*89;
  "qbits" = 3;
  chute =2;
  q = 2^20;
Tempo de execução:
  Cuda: 0.15 segundos Nao acertou, chute talvez
  C: 0.036610 segundos
  Python: 0.3879 segundos
  
Quinto teste
  N = 47*31;
  "qbits" = 3;
  chute =2;
  q = 2^20;
Tempo de execução:
  Cuda: 0.159 segundos
  C: 0.961187 segundos
  Python: 0.3933 segundos

Sexto teste
  N = 47*31;
  "qbits" = 3;
  chute =2;
  q = 2^24;
Tempo de execução:
  Cuda:  3.741 segundos
  C: 0.045941 segundos Nao acertou, chute talvez
  Python: 6.1758 segundos
  
Setimo teste
  N = 59*41;
  "qbits" = 3;
  chute =2;
  q = 2^20;
Tempo de execução:
  Cuda:  0.16 segundos
  C: 0.056414 segundos Nao acertou, chute talvez
  Python: 0.3919 segundos
  
Oitavo teste
  N = 59*41;
  "qbits" = 3;
  chute =2;
  q = 2^24;
Tempo de execução:
  Cuda:  0.996 segundos Nao acertou, chute talvez
  C: 0.936393 segundos
  Python: 5.9983 segundos
  
Nono teste
  N = 23*41;
  "qbits" = 3;
  chute =2;
  q = 2^24;
Tempo de execução:
  Cuda:  2.054 segundos Nao acertou, chute talvez
  C: 0.940327 segundos Nao acertou, chute talvez
  Python: 6.1717 segundos
  
Décimo teste
  N = 23*41;
  "qbits" = 3;
  chute =2;
  q = 2^20;
Tempo de execução:
  Cuda:  0.173 segundos
  C: 0.040281 segundos Nao acertou, chute talvez
  Python: 0.4062 segundos

  Valores testado para permutação de 23*29, 23*31, 23*37, 23*41, 23*47, 23*53, 23*59, 23*61, 23*67, 23*71, 23*73, 23*79, 23*83, 23*89, 23*97
  Valores testado para permutação de 31*23,31*23, 31*37, 31*41, 31*47, 31*53, 31*59, 31*61, 31*67, 31*71, 31*73, 31*79, 31*83, 31*89, 31*97
  Valores testado para permutação de 29*23, 29*31, 29*37, 29*41, 29*47, 29*53, 29*59, 29*61, 29*67, 29*71, 29*73, 29*79, 29*83, 29*89, 29*97

  todos os 3 testes nao obteve todos resultados, tanto em python, quanto em C, quando em CUDA.
  Falha no estouro no mmc erra toda a conta para achar o valor correto.
  o proprio python ele estoura, porém ha tratamento de erro que é manter o ultimo numero que nao foi estourado. Não sei o motivo de dar certo ainda.
  passa o algoritmo de frações continuadas para long int e obteve resultado melhor que o python. porém na hora de aplicar existe erro ainda. 
  função fat, as vezes so não funciona. Usei printf para poder achar o erro, ele entra nos laçoes e condinções porem ele nao executa outros comandos. Provavel ser memory leak.
  Preciso terminar esse código, deixa-lo o mais otimizado possivel. Em CUDA não consegui achar uma solução melhor a não ser que desparalelize algumas funções, tirando isso, tem o mesmo "defeito" que C, porém defeito ja descoberto e tendo tentativas de soluções.

  Para poder compilar seu código, é necessário ter a biblioteca fftw3 instalado e o comando do terminal ficará assim:
  gcc teste.c -o teste -lfftw3 -lgmp -lm

  Saidas diferentes para o calculo da probabilidade em c e python. alem disso resultado das medidas muitos diferentes, porém por ser número aleatório. Dentro do limite

  
![image](https://github.com/raianffc/Simulacao-Quantica-usando-GPU/assets/54862169/b51d44e9-316d-4dac-9375-984e149eb2dd)
Nesse print a probababilidade saiu igual!
Result: 
C
[436908 961195 576718 506813 716528 850512 990323 716528 256320 413605 168937 984496 1 145636 891291 ]
Python
[343700, 774782, 570891, 367001, 93206, 1042750, 425256, 1, 757304, 1042750, 343700, 1042750, 297096, 75730, 413604]
