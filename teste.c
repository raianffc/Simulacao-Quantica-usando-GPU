#include <stdio.h>
#include <complex.h>
#include <math.h>

#define PI 3.14159265358979323846

static void swap(double complex *a, double complex *b) {
    double complex temp = *a;
    *a = *b;
    *b = temp;
}

void fft(double complex *data, int size, int inverse) {
    int bits, i, j, k;
    double complex w, w_m, temp;

    j = 0;
    for (i = 0; i < size - 1; i++) {
        if (i < j)
            swap(&data[i], &data[j]);
        k = size / 2;
        while (k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }

    for (bits = 1; bits < size; bits *= 2) {
        double angle = (2.0 * PI / bits) * (inverse ? 1 : -1);
        w_m = cos(angle) + sin(angle) * I;
        w = 1.0 + 0.0 * I;
        for (j = 0; j < bits / 2; j++) {
            for (i = j; i < size; i += bits) {
                k = i + bits / 2;
                temp = w * data[k];
                data[k] = data[i] - temp;
                data[i] += temp;
            }
            w *= w_m;
        }
    }
}

void ifft(double complex *data, int size) {
    fft(data, size, 1);
    for (int i = 0; i < size; ++i) {
        data[i] /= size;
    }
}

int main() {
    int size = 2;
    double complex Z[] = {0 + 0*I, 1 + 0*I};

    printf("Dados originais:\n");
    for (int i = 0; i < size; ++i) {
        printf("%f + %fi\n", creal(Z[i]), cimag(Z[i]));
    }

    fft(Z, size, 0); // Calcula a FFT

    printf("\nResultado da FFT:\n");
    for (int i = 0; i < size; ++i) {
        printf("%f + %fi\n", creal(Z[i]), cimag(Z[i]));
    }

    ifft(Z, size); // Calcula a IFFT

    printf("\nResultado da IFFT:\n");
    for (int i = 0; i < size; ++i) {
        printf("%f + %fi\n", creal(Z[i]), cimag(Z[i]));
    }

    return 0;
}