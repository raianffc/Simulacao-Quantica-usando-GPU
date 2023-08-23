#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <curand_kernel.h>

#define PI 3.14159265358979323846

__global__ void generateRandom(float *randomValues, unsigned int seed) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    curandState state;
    curand_init(seed, tid, 0, &state);

    randomValues[tid] = curand_uniform(&state);
}

int main() {
    int numThreads = 1;
    int numBlocks = 1;
    int totalThreads = numThreads * numBlocks;

    float *d_randomValues;
    float *h_randomValues = (float *)malloc(totalThreads * sizeof(float));

    cudaMalloc((void **)&d_randomValues, totalThreads * sizeof(float));

    generateRandom<<<numBlocks, numThreads>>>(d_randomValues, time(0));
    
    cudaMemcpy(h_randomValues, d_randomValues, totalThreads * sizeof(float), cudaMemcpyDeviceToHost);

    for (int i = 0; i < totalThreads; i++) {
        printf("Random Value %d: %f\n", i, h_randomValues[i]);
    }

    cudaFree(d_randomValues);
    free(h_randomValues);

    return 0;
}

