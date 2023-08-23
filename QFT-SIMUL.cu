#include <stdio.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>
#include <stdlib.h>
#include <math.h>
#include <curand_kernel.h>
#include<stdbool.h>
#include<time.h>
#include<complex.h>

int mdc(int num1, int num2) {
    int resto;
    do {
        resto = num1 % num2;
        num1 = num2;
        num2 = resto;
    } while (resto != 0);
    return num1;
}

int mmc(int num1, int num2) {
    int a;
    if(num2==0) return num1;
    a = mdc(num1,num2);
    return (num1 * num2) / a;
    
}
float buscabin(float *Soma, float *P, double m, int tamSoma) {
    int n = tamSoma;
    if (n == 0) {
        return 0;
    } else if (n == 1) {
        return P[0];
    } else if (n == 2) {
        if (m <= Soma[0]) {
            return P[0];
        } else {
            return P[1];
        }
    } else {
        int meio = n / 2;
        if (m == Soma[meio]) {
            return P[meio];
        } else if (m < Soma[meio]) {
            return buscabin(Soma, P, m, meio);
        } else {
            return buscabin(&Soma[meio], &P[meio], m, n - meio);
        }
    }
}
double* FracCont(double x, double q, double N, int *tamL) {
    int tam = 1;
    double x_inic = x;
    double *L;
    L = (double*)malloc(tam * sizeof(double));
    if (L == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    if (x == 1) {
        tam = 2;
        L = (double*)realloc(L, tam * sizeof(double));
        if (L == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
        L[0] = 1;
        L[1] = 0;
        x = 0;
    } else {
        L[0] = -1;
    }
    x = x / q;
    int i = 0;
    double max = log(x_inic) / log(1.6);
    do {
        int c = (int) x;
        if (tam != 1) {
            L = (double*)realloc(L, (tam + 1) * sizeof(double));
            if (L == NULL) {
                printf("Erro na alocacao de memoria.");
                exit(1);
            }
            L[tam] = (double) c;
        } else {
            L[0] = (double) c;
        }
        x = x - c;
        if (x >= 0) {
            x = 1 / x;
            if (x > x_inic || x > N) {
                x = 0;
            }
        } else {
            x = 0;
        }
        i++;
        tam++;
    } while (x > 0 && i < max);
    *tamL = tam - 1;
    return L;
}

double** Frac(double *L, int *tamL) {
    int tamF = *tamL;
    if (*tamL == 0) {
        *tamL = 1;
        L = (double*)realloc(L, (*tamL) * sizeof(double));
        if (L == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
        L[0] = 1;
    }
    double **F = (double**)malloc(2 * sizeof(double *));
    if (F == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for (int i = 0; i < 2; i++) {
        F[i] = (double*)malloc((tamF) * sizeof(double));
        if (F[i] == NULL) {
            printf("\nErro na alocacao de memoria.\n");
            exit(1);
        }
    }
    for(int i = 0;i<2; i++){
        for(int j = 0; j<(tamF);j++){
            F[i][j]=0.0;
        }
    }
    int num = 1;
    int den = 1;
    for ( int i=0; i<((tamF)-1);i++){
        num = 1;
        den = (int)L[(tamF)-i-1];
        for (int j=((tamF)-i-2); j>0;j--){
            int temp=num;
            num = den;
            den = L[j]*den+temp;
        }
        F[0][i+1]=(double)num;
        F[1][i+1]=(double)den;
    }
    F[0][0]=(tamF)-1;

    int total = (int)F[1][2];
    for(int j=3;j<tamF;j++){
        total = mmc(total,(int)F[1][j]);
    }
    F[1][0]=(double)total;
    if(F[1][0]==0 && *tamL==1){
        printf("\nFracao inexistente: %.0f \n", L[0]);
    }
    free(L);
    return F;
}

__global__ void calculateZ(cufftComplex *d_Y, int q, int r) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < q) {
        d_Y[idx] = make_cuComplex(0.0f, 0.0f);
    }

    __syncthreads();

    int j = 1;
    while (j <= q) {
        if (idx == j) {
            d_Y[j] = make_cuComplex(1.0f, 0.0f);
        }
        j += r;
        __syncthreads();
    }
}

__global__ void Quadrado_conjugado(double *d_Z, cufftComplex *d_Y, int q) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < q) {
        d_Z[idx] = cuCabsf(d_Y[idx]) * cuCabsf(d_Y[idx]);
    }
}

__global__ void normalizeZ(double *d_Z, double sum_Z, int q) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < q) {
        d_Z[idx] /= sum_Z;
    }
}
__global__ void sumVector(float *input, int size, float *result) {
    extern __shared__ float sdata[];
    
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    sdata[tid] = (i < size) ? input[i] : 0.0f;
    
    __syncthreads();
    
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) {
        result[blockIdx.x] = sdata[0];
    }
}


double *Prepara(double N, double x, double *r, double q) {
    int tamN = (int)log2(N);
    double q1 = 1 << (2 * tamN);  
    printf("Valor ideal para q: %.0f\n", q1);

    if (q < N) {
        q = 1 << (tamN + 4);
    }

    if (*r == 0) {
        int s = x;
        int i = 1;
        while (s > 1) {
            s = (int)(s * x) % ((int)N);
            i++;
        }
        *r = i;
        printf("Ordem r nao informada. Ordem r calculada: %.0f\n", *r);
    } else {
        printf("Ordem r informada: %f\n", *r);
    }

    printf("Criando Z...\n");
    double *Z;
    cudaMallocHost((void **)&Z, q * sizeof(double));

    // Aloca memória na GPU para Y
    cufftComplex *d_Y;
    cudaMalloc((void **)&d_Y, q * sizeof(cufftComplex));

    // Define o número de threads por bloco e calcula o número de blocos
    int threadsPerBlock = 256; 
    int numBlocks = (q + threadsPerBlock - 1) / threadsPerBlock;

    // Preenche Y
    calculateZ<<<numBlocks, threadsPerBlock>>>(d_Y, q, (int)(*r));
    cudaDeviceSynchronize();

    // Calcula FFT
    printf("Calculando FFT...\n");
    cufftHandle plan;
    cufftPlan1d(&plan, q, CUFFT_C2C, 1);
    cufftExecC2C(plan, d_Y, d_Y, CUFFT_FORWARD);
    cufftDestroy(plan);

    // Calcula probabilidades
    printf("Calculando probabilidades...\n");
    Quadrado_conjugado<<<numBlocks, threadsPerBlock>>>(Z, d_Y, q);
    cudaDeviceSynchronize();

    // Calcula a soma das probabilidades
    double sum_Z = 0;
    for (int i = 0; i < q; i++) {
        sum_Z += Z[i];
    }
    double temp=0;
    for (int i = 0; i < q; i++) {
        Z[i] = Z[i]/sum_Z;
        temp += Z[i];
    }
    sum_Z = temp;
    printf("Soma das probabilidades: %.20f\nCriando Soma com probabilidade acumulada...\n", sum_Z);
    


    cudaFree(d_Y);
    return Z;
}
/*__global__ void calculateSomaPParallel(float *P, float *Soma, double *Z, double k, double q, int r, curandState *state) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < r) {
        double pos = idx * k;

        P[2 * idx] = (float)pos;
        P[2 * idx + 1] = pos + 1;

        double total = 0.0;
        for (int i = 0; i <= pos; i++) {
            total += Z[i];
        }

        // Gerar um número aleatório entre 0 e 100 usando o cuRAND
        float random_value = curand_uniform(&state[idx]) * 101;

        Soma[2 * idx] = total;
        Soma[2 * idx + 1] = total + Z[(int)pos + 1] * random_value;
    }
}

float *Soma_P_Parallel(double r, double q, float *P, float *Soma, double *Z) {
    double k = (q / r);

    // Aloca estados do cuRAND
    int numThreadsPerBlock = 256;
    int numBlocks = (r + numThreadsPerBlock - 1) / numThreadsPerBlock;
    curandState *devStates;
    cudaMalloc((void **)&devStates, numBlocks * numThreadsPerBlock * sizeof(curandState));
    setupCurand<<<numBlocks, numThreadsPerBlock>>>(devStates, time(0));
    cudaDeviceSynchronize();

    // Configuração dos blocos e threads para o kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (r + threadsPerBlock - 1) / threadsPerBlock;

    // Chama o kernel para calcular Soma_P paralelamente
    calculateSomaPParallel<<<blocksPerGrid, threadsPerBlock>>>(P, Soma, Z, k, q, r, devStates);
    cudaDeviceSynchronize();

    // Libera memória alocada para os estados do cuRAND
    cudaFree(devStates);

    return Soma;
}*/ //Nao sabendo calcular Probabilidade Acumulada (dos picos)
__global__ void generateRandom(float *randomValues, unsigned int seed) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    curandState state;
    curand_init(seed, tid, 0, &state);

    randomValues[tid] = curand_uniform(&state);
}
float *Soma_P (double r, double q, float *P, float *Soma, double *Z){
    double k = (q / r);
    printf("Calculando Somas...\n");
    double total = 0;

    for (int i = 0; i <r; i++) {
        double pos = (i * k);
        P[2 * i] = (float)pos;
        total += Z[(int)pos];
        Soma[2*i]= total;
        P[2*i+1]= pos+1;
        total = total + Z[((int)pos)+1];
        Soma[2*i+1] = total;
    }
    int numThreads = 1;
    int numBlocks = 1;
    int totalThreads = numThreads * numBlocks;

    float *d_randomValues;
    float *h_randomValues = (float *)malloc(totalThreads * sizeof(float));

    cudaMalloc((void **)&d_randomValues, totalThreads * sizeof(float));

    generateRandom<<<numBlocks, numThreads>>>(d_randomValues, time(0));
    
    cudaMemcpy(h_randomValues, d_randomValues, totalThreads * sizeof(float), cudaMemcpyDeviceToHost);

    P[2*(((int)r)-1)] = -1 * (int)(h_randomValues[0] * 101 * q);
    Soma[2*(((int)r)-1)] = 1;

    printf("Probabilidade Acumulada (dos picos): %f\n", total);
    
    cudaFree(d_randomValues);
    free(h_randomValues);

    return Soma;
}
double *Simula(float *Soma,float *P, int tamResult, int tamSoma){
    srand(time(NULL));
    printf("\nSimula medicao:\n");
    double *result;
    result=(double*)malloc(tamResult*sizeof(double));
    if (result == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for(int i = 0; i<tamResult; i++) result[i]=0;
    int numThreads = 256; 
    int numBlocks = (tamResult + numThreads - 1) / numThreads;
    int totalThreads = numThreads * numBlocks;

    float *d_randomValues;
    float *h_randomValues = (float *)malloc(totalThreads * sizeof(float));

    cudaMalloc((void **)&d_randomValues, totalThreads * sizeof(float));

    generateRandom<<<numBlocks, numThreads>>>(d_randomValues, time(0));
    
    cudaMemcpy(h_randomValues, d_randomValues, totalThreads * sizeof(float), cudaMemcpyDeviceToHost);

    for(int i = 0; i<tamResult; i++){
    // Normaliza o número para estar entre 0 e 1, incluindo valores decimais
        double m = (double)h_randomValues[i];
        result[i] = (double)buscabin(Soma,P,m, tamSoma);
        if (result[i]==0){
            result[i] = 1;
        }
    }

    cudaFree(d_randomValues);
    free(h_randomValues);

    return result;
}
float **EstimaOrdem(double r,double *result,double q, double N, int n){
    float **R;
    int taml = 1;
    R=(float**)malloc(n*sizeof(float*));
    if (R == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for(int i=0;i<n;i++){
        R[i]=(float*)malloc(3*sizeof(float));
        if (R[i] == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
    }
    printf("\nTenta estimar a ordem r=ord(x,N) ou multiplo ou divisor dela para extrair os fatores de N\n");
    for (int i=0;i<n;i++){
        for(double j =-1; j<2; j++){
            printf(". ");
            double *l = FracCont(result[i]+j,q,N, &taml);
            double **t = Frac(l, &taml);
            int k = (float)j;
            R[i][k+1]=(float)t[1][0];
            free(t);
        }
    }
    return R;
}
float **EstimaFator(double N, double x,float **R, int tam){
    printf("\nProcura um multiplo da ordem ou um divisor que distiga um fator nao trivial.\n");
    float **Sucesso;
    int potTotal=1;
    int pot;
    Sucesso=(float**)malloc(tam*sizeof(float*));
    if (Sucesso == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for(int i = 0; i<tam; i++){
        Sucesso[i]=(float*)malloc(3*sizeof(float));
        if (Sucesso[i] == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
    }
    for(int i = 0; i<tam; i++){
        float sucesso[3]={0,0,0};
        for(int j=0; j<3; j++){
            double x1=(double)x;
            double R1=(double)R[i][j];
            pot=(((int)(pow(x1,R1)))%((int)N));
             if (pot==1){
                sucesso[j]=1;
                printf("\nmultiplo de r\n");
            }
            else{
                int d=mdc(pot-1,N);
                if (d > 1){
                    sucesso[j]=2;
                    printf("\nfator: %d\n",d);
                }
            }
            potTotal = potTotal*pot;
            if(potTotal%((int)N)==1){
                sucesso[j]=sucesso[j]+4;
                printf("\nfatores do múltiplo da ordem %d %d\n",pot,(potTotal/pot));
            }
            else{
                int d=mdc((potTotal-1),N);
                if (d > 1){
                    sucesso[j]=sucesso[j]+8;
                }
            }
            potTotal = potTotal%((int)N);
            //printf("potTotal: %d", potTotal);
            Sucesso[i][j]=sucesso[j];
        }

    }

    return Sucesso;
}
int existe(int valor, int *array, int tamanho) {
    for (int i = 0; i < tamanho; i++) {
        if (array[i] == valor) {
            return 1;  // O valor já existe
        }
    }
    return 0;  // O valor não existe
}
int removerDuplicatas(int *array, int tamanho) {
    if (tamanho <= 1) {
        return tamanho;  // Não há duplicatas para remover
    }

    int novoTamanho = 1;  // Tamanho do novo array sem duplicatas
    for (int i = 1; i < tamanho; i++) {
        if (!existe(array[i], array, novoTamanho)) {
            array[novoTamanho] = array[i];  // Adiciona o elemento único
            novoTamanho++;
        }
    }

    return novoTamanho;
}
int* Fatores(double N, double x, float **R, float **S, int num_s, int *num_fatores) {
    int *fat;
    int *aux;
    fat = (int *)malloc((num_s * num_s * 2) * sizeof(int));
    if (fat == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    int count = 0;

    for (int i = 0; i < num_s; i++) {
        for (int j = 0; j < 3; j++) {
            if (S[i][j] == 1) {
                if ((int)R[i][j] % 2 == 0) {
                    int f = mdc((int)(((int)pow(x, (int)(R[i][j] / 2))%(int)N) - 1), (int)N);
                    fat[count++] = f;
                }
                if ((int)R[i][j] % 3 == 0) {
                    int f = mdc((int)(((int)pow(x, (int)(R[i][j] / 3))%(int)N) - 1), (int)N);
                    fat[count++] = f;
                }
            } else if (S[i][j] == 2) {
                int f = mdc((int)pow(x, (int)R[i][j]) - 1, (int)N);
                fat[count++] = f;
                fat[count++] = (int)N / f;
            }
        }
    }
   /*int k=1, j=0;
    for(int i=0; i<count;i++){
        if(fat[i]!=N){
            aux = (int*)realloc(aux,k*sizeof(int));
            aux[j] =fat[i];
		    k++;
            j++;
        }
    }*/
    //printf("\ncontador k: %d \n",k);
    /*if(k>1){
        printf("Aux: [");
        for(int i=0; i<k;i++){
            printf("%d ", aux[i]);
        }
        printf("]\n");
    }*/
    *num_fatores = removerDuplicatas(fat, count);
    return fat;
}
int main(){

    //double time_spent = 0.0;
    //clock_t begin = clock();
    double p1 = 31;
    double p2 = 29;
    double N  = p1 * p2; //N nao precisa ser semi-primo
    double x  = 2;
    double r  = 0;
    double q  = (int)pow(2, 24);//2**20
    int n  = 15; // quantidade de valores medidos 
    float *Soma;
    float *P;
    int tamSoma_P;
    float **R;
    float **S;
    int *fat;
    int tamFat;
    double *Z;
    int threadsPerBlock = 256;
    int numBlocks = ( (int)r+ threadsPerBlock - 1) / threadsPerBlock;

    Z = Prepara(N, x, &r, q);

    int r_int = (int)r;

    P = (float *)malloc((2 * (r + 1)) * sizeof(float));
    Soma = (float *)malloc((2 * r + 1) * sizeof(float));
    tamSoma_P=(2*(r+1));
    if (P == NULL || Soma == NULL || Z == NULL) {
        printf("Erro na alocação de memória.");
        exit(1);
    }

    Soma = Soma_P(r, q, P, Soma, Z);

    printf("%.0f", r);

    double *result;
    result = Simula(Soma, P, n, tamSoma_P);
    printf("\nq= %.0f\n",q);
    printf("\nResultados medidos na rotina quantica do QOFA:\n");
    printf("[");
    for(int i=0; i<n;i++){
        printf("%.0f ", result[i]);
    }
    printf("]\n");

    R= EstimaOrdem(r, result, q, N, n);

    S = EstimaFator(N, x, R, n);
    fat = Fatores(N, x, R, S, n, &tamFat);
    printf("\nFatores: \n");
    printf("[");
    for(int i=0; i<tamFat;i++){
        printf("%d ",fat[i]);
    }
    printf("]\n");
    free(fat);
    
    cudaFreeHost(Z);
    free(P);
    free(Soma);
    //clock_t end = clock();
    //time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
 
    //printf("\nTempo de execucao: %f segundos\n", time_spent);
    return 0;
}