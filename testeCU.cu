#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include<stdbool.h>
#include<math.h>
#include<time.h>
#include<complex.h>
#include <fftw3.h>
#include <cuComplex.h>
#include <cufft.h>

__global__ void calculateProbabilities(cuDoubleComplex *Y, double *Z, int q) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < q) {
        Z[idx] = cuCabs(Y[idx]) * cuCabs(Y[idx]);
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
    double *Z = (double *)malloc(q * sizeof(double));
    if (Z == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }

    cuDoubleComplex *Y;
    cudaMalloc((void **)&Y, q * sizeof(cuDoubleComplex));
    cudaMemset(Y, make_cuDoubleComplex(0.0, 0.0), q * sizeof(cuDoubleComplex));

    calculateProbabilities<<<(q + 255) / 256, 256>>>(Y, Z, q);

    cudaFree(Y);

    double sum_Z = 0;
    for (int i = 0; i < q; i++) {
        sum_Z += Z[i];
    }

    double temp = 0;
    for (int i = 0; i < q; i++) {
        Z[i] = Z[i] / sum_Z;
        temp += Z[i];
    }
    sum_Z = temp;
    printf("Soma das probabilidades: %.20f\nCriando Soma com probabilidade acumulada...\n", sum_Z);

    return Z;
}
__global__ void calculateSums(float *P, float *Soma, double *Z, int r, double q) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < r) {
        double k = (q / r);
        double pos = (idx * k);

        P[2 * idx] = (float)pos;
        Soma[2 * idx] = Z[(int)pos];

        if (idx > 0) {
            Soma[2 * idx] += Soma[2 * (idx - 1)];
        }

        P[2 * idx + 1] = pos + 1;

        if (idx < r - 1) {
            Soma[2 * idx + 1] = Soma[2 * idx];
        } else {
            P[2 * idx + 1] = -1 * (int)(curand_uniform(&state) * 101 * q);
            Soma[2 * idx + 1] = 1;
        }
    }
}

__global__ void Soma_P_CUDA(double r, double q, float *P, float *Soma, double *Z, int n) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < n) {
        double k = (q / r);
        double pos = tid * k;
        P[2 * tid] = (float)pos;
        atomicAdd(&Soma[2 * tid], Z[(int)pos]);
        P[2 * tid + 1] = pos + 1;
        atomicAdd(&Soma[2 * tid + 1], Z[((int)pos) + 1]);

        if (tid == ((int)r) - 1) {
            P[2 * tid] = -1 * (int)((rand() % 101) * q);
            Soma[2 * tid] = 1;
        }
    }
}

float *Soma_P_CUDA_Wrapper(double r, double q, float *P, float *Soma, double *Z, int n) {
    float *d_P, *d_Soma;
    double *d_Z;

    cudaMalloc((void **)&d_P, 2 * n * sizeof(float));
    cudaMalloc((void **)&d_Soma, 2 * n * sizeof(float));
    cudaMalloc((void **)&d_Z, q * sizeof(double));

    cudaMemcpy(d_P, P, 2 * n * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Soma, Soma, 2 * n * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Z, Z, q * sizeof(double), cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;
    int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

    Soma_P_CUDA<<<blocksPerGrid, threadsPerBlock>>>(r, q, d_P, d_Soma, d_Z, n);

    cudaMemcpy(P, d_P, 2 * n * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(Soma, d_Soma, 2 * n * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(d_P);
    cudaFree(d_Soma);
    cudaFree(d_Z);

    return Soma;
}
int main(){
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
    
    Z = Prepara(N, x, &r, q);
    cudaDeviceSynchronize();
    tamSoma_P=(2*(r+1));
    P = (float*)malloc((2*(r+1))*sizeof(float));
    Soma = (float*)malloc((2*r+1)*sizeof(float));
    if (P == NULL || Soma == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for(int i=0;i<(2*(r+1));i++){
        P[i]=0;
        Soma[i]=0;
    }
    Soma = Soma_P(r, q, P, Soma, Z);
    
    
    printf("%.0f", r);
    
    return 0;
}