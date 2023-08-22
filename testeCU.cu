#include <stdio.h>
#include <cuda_runtime.h>

__global__ void vectorAdd(const double *A, const double *B, double *C, int size) {
    printf("\naqui2\n");
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < size) {
        C[tid] = A[tid] + B[tid];
    }
}

int main() {
    int size = 5; // Tamanho dos vetores
    double *h_A, *h_B, *h_C; // Vetores na CPU
    double *d_A, *d_B, *d_C; // Vetores na GPU

    // Aloca memória para os vetores na CPU
    h_A = (double *)malloc(size * sizeof(double));
    h_B = (double *)malloc(size * sizeof(double));
    h_C = (double *)malloc(size * sizeof(double));

    // Inicializa os vetores na CPU
    for (int i = 0; i < size; i++) {
        h_A[i] = i;
        h_B[i] = 2 * i;
    }

    // Aloca memória para os vetores na GPU
    cudaMalloc((void **)&d_A, size * sizeof(double));
    cudaMalloc((void **)&d_B, size * sizeof(double));
    cudaMalloc((void **)&d_C, size * sizeof(double));

    // Copia os vetores da CPU para a GPU
    cudaMemcpy(d_A, h_A, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, size * sizeof(double), cudaMemcpyHostToDevice);

    // Define o número de threads por bloco
    int threadsPerBlock = 256;
    int numBlocks = (size + threadsPerBlock - 1) / threadsPerBlock;

    // Executa o kernel de soma dos vetores
    printf("\naqui1\n");
    vectorAdd<<<numBlocks, threadsPerBlock>>>(d_A, d_B, d_C, size);
    cudaDeviceSynchronize();
    printf("\naqui3\n");
    // Copia o vetor resultado da GPU para a CPU
    cudaMemcpy(h_C, d_C, size * sizeof(double), cudaMemcpyDeviceToHost);

    // Exibe o resultado
    printf("Vetor A: ");
    for (int i = 0; i < size; i++) {
        printf("%.1f ", h_A[i]);
    }
    printf("\nVetor B: ");
    for (int i = 0; i < size; i++) {
        printf("%.1f ", h_B[i]);
    }
    printf("\nResultado da soma: ");
    for (int i = 0; i < size; i++) {
        printf("%.1f ", h_C[i]);
    }
    printf("\n");

    // Libera a memória alocada
    free(h_A);
    free(h_B);
    free(h_C);
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    cudaError_t cudaError = cudaGetLastError();
    if (cudaError != cudaSuccess) {
        printf("Erro durante a execução do kernel: %s\n", cudaGetErrorString(cudaError));
        return 1; // Encerra o programa com um código de erro
    }

    return 0;
}
