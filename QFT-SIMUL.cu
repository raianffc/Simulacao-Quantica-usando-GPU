#include <stdio.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>
#include <stdlib.h>
#include <math.h>

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

__global__ void calculateProbabilities(double *d_Z, cufftComplex *d_Y, int q) {
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
__global__ void reduceSum(double *d_input, double *d_output, int N) {
    extern __shared__ double shared_data[]; // Memória compartilhada

    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + tid;

    if (idx < N) {
        shared_data[tid] = d_input[idx];
    } else {
        shared_data[tid] = 0.0;
    }

    __syncthreads();

    // Realiza a redução na memória compartilhada
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            shared_data[tid] += shared_data[tid + s];
        }
        __syncthreads();
    }

    // O thread 0 de cada bloco escreve o resultado final no vetor de saída
    if (tid == 0) {
        d_output[blockIdx.x] = shared_data[0];
    }
}

__global__ void normalizeZ(double *d_Z, double sum_Z, int q) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < q) {
        d_Z[idx] /= sum_Z;
    }
}
__global__ void parallelSoma_P(float *d_P, float *d_Soma, double *d_Z, int r, int q) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < r) {
        double k = (double)q / r;
        double pos = idx * k;

        d_P[2 * idx] = (float)pos;
        double total = 0.0;
        for (int i = 0; i <= idx; i++) {
            total += d_Z[(int)(i * k)];
        }
        d_Soma[2 * idx] = total;
        d_P[2 * idx + 1] = (float)pos + 1;
        d_Soma[2 * idx + 1] = total + d_Z[(int)((idx + 1) * k)];
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
        printf("Ordem r não informada. Ordem r calculada: %.0f\n", *r);
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
    calculateProbabilities<<<numBlocks, threadsPerBlock>>>(Z, d_Y, q);
    cudaDeviceSynchronize();

    // Calcula a soma das probabilidades
    double *d_Z;
    double sum_Z;    
    cudaMalloc((void **)&d_Z, q * sizeof(double));

    // Copia os dados da CPU para a GPU
    cudaMemcpy(d_Z, Z, q * sizeof(double), cudaMemcpyHostToDevice);
    int sharedMemorySize = threadsPerBlock * sizeof(double);

    // Executa o kernel de redução para calcular a soma
    reduceSum<<<numBlocks, threadsPerBlock, sharedMemorySize>>>(d_Z, &sum_Z, q);
    cudaDeviceSynchronize();

    // Executa o kernel para normalizar Z
    normalizeZ<<<numBlocks, threadsPerBlock>>>(d_Z, sum_Z, q);
    cudaDeviceSynchronize();
    //calcula a probabilide entre 0 e 1
    reduceSum<<<numBlocks, threadsPerBlock, sharedMemorySize>>>(d_Z, &sum_Z, q);
    cudaDeviceSynchronize();
    // Copia o resultado de volta para a CPU
    cudaMemcpy(Z, d_Z, q * sizeof(double), cudaMemcpyDeviceToHost);
    printf("Soma das probabilidades: %.20f\nCriando Soma com probabilidade acumulada...\n", sum_Z);


    cudaFree(d_Y);
    cudaFree(d_Z);
    return Z;
}

float *ParallelSoma_P(double r, double q, float *P, float *Soma, double *Z) {
    int r_int = (int)r;
    int threadsPerBlock = 256;
    int numBlocks = (r_int + threadsPerBlock - 1) / threadsPerBlock;

    float *d_P, *d_Soma;
    double *d_Z;

    cudaMalloc((void **)&d_P, (2 * r_int) * sizeof(float));
    cudaMalloc((void **)&d_Soma, (2 * r_int) * sizeof(float));
    cudaMalloc((void **)&d_Z, q * sizeof(double));
    
    cudaMemcpy(d_Z, Z, q * sizeof(double), cudaMemcpyHostToDevice);

    parallelSoma_P<<<numBlocks, threadsPerBlock>>>(d_P, d_Soma, d_Z, r_int, (int)q);

    cudaMemcpy(P, d_P, (2 * r_int) * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(Soma, d_Soma, (2 * r_int) * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(d_P);
    cudaFree(d_Soma);
    cudaFree(d_Z);

    return Soma;
}


int main() {
    double p1 = 31;
    double p2 = 29;
    double N = p1 * p2;
    double x = 2;
    double r = 0;
    double q = pow(2, 20);
    double *Z;
    float *P;
    float *Soma;
    int threadsPerBlock = 256;
    int numBlocks = ( (int)r+ threadsPerBlock - 1) / threadsPerBlock;

    Z = Prepara(N, x, &r, q);

    int r_int = (int)r;
    P = (float *)malloc(2 * r_int * sizeof(float));
    Soma = (float *)malloc(2 * r_int * sizeof(float));
    ParallelSoma_P(r, q, P, Soma, Z);    
    
    
    cudaFreeHost(Z);
    cudaDeviceSynchronize();
    return 0;
}