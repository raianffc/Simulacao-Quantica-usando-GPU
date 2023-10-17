#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include<stdbool.h>
#include<math.h>
#include<time.h>
#include<complex.h>
#include <fftw3.h>
#include <unistd.h> 
#define PI 3.14159265358979323846
double gcd(double a, double b) {
    while (b != 0) {
        double temp = b;
        b = fmod(a, b);
        a = temp;
    }
    return a;
}

// Função para calcular o LCM (Mínimo Múltiplo Comum)
double lcm(double a, double b) {
    return (a * b) / gcd(a, b);
}
long int mdc(long int num1, long int num2) {
    long int resto;
    do {
        resto = num1 % num2;
        num1 = num2;
        num2 = resto;
    } while (resto != 0);
    return num1;
}

long int mmc(long int num1, long int num2) {
    long int resto, a;
    if(num2==0) return num1;
    a = mdc(num1,num2);
    return (num1 * num2) / a;
    
}
long int* FracCont(double x, double q, double N, int* length) {
    x = fabs(x);
    double x_inic = x;
    long int* L;
    if (x == 1.0) {
        L = (long int*)malloc(2 * sizeof(long int));
        if (L == NULL) {
            perror("Erro ao alocar memória");
            exit(1);
        }
        L[0] = 1;
        L[1] = 0;
        x = 0;
        *length = 2;
        return L;
    } else {
        long int* L = NULL;
        *length = 0;
    }
    x = x / q;
    int i = 0;
    double max = log(x_inic) / log(1.6);
    while ((x > 0.0 && i < max) && (x_inic !=1)) {
        int c = (int)x;
        (*length)++;
        long int* temp = (long int*)realloc(L, (*length) * sizeof(long int));
        if (temp == NULL) {
            perror("Erro ao alocar memória");
            exit(1);
        }
        L = temp;
        L[(*length) - 1] = (long int)c;
        x -= c;
        if (x >= 0.0) {
            x = 1.0 / x;
            if (x > x_inic || x>N) {
                x = 0.0;
            }
        } else {
            x = 0.0;
        }
        i++;
    }
    return L;
}
long int** Frac(long int* L, int m) {
    long int** F = (long int**)malloc(2 * sizeof(long int*));
    if (F == NULL) {
        perror("Erro ao alocar memória");
        exit(1);
    }
    for (int i = 0; i < 2; i++) {
        F[i] = (long int*)malloc(m * sizeof(long int));
        if (F[i] == NULL) {
            perror("Erro ao alocar memória");
            exit(1);
        }
    }
    long int num = 1.0;
    long int den = 1.0;
    if (m == 0) {
        m = 1;
        long int* temp = (long int*)malloc(sizeof(long int));
        if (temp == NULL) {
            perror("Erro ao alocar memória");
            exit(1);
        }
        temp[0] = 1;
        free(L);
        L = temp;
    }
    for (int i = 0; i < m - 1; i++) {
        num = 1.0;
        den = L[m - i - 1];
        for (int j = m - i - 2; j > 0; j--) {
            long int temp_num = num;
            num = den;
            den = L[j] * den + temp_num;
        }
        F[0][i + 1] = num;
        F[1][i + 1] = den;
    }
    F[0][0] = m - 1;
    if (F[1][2] != 0.0) {
        if (m > 1) {
            // Você precisará implementar a função "lcm" ou usar uma implementação existente.
            long int total = F[1][2];
           // long int maior=0;
            for(int j=3;j<m;j++){
            total = mmc(total,F[1][j]);
            //if(maior<total) maior=total;
            }
            if(total<0) total=total*-1;
            F[1][0]=total;
            printf("\nFzin: %ld\n",F[1][0]);

        }
    }
    if (F[1][1] ==0) {
        printf("Fração inexistente:");
    }

    return F;
}
long int **EstimaOrdem(double r,double *result,double q, double N, int n){
    int mult = 0;
    long int **R;
    int taml;
    R=(long int**)malloc(n*sizeof(long int*));
    if (R == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for(int i=0;i<n;i++){
        R[i]=(long int*)malloc(3*sizeof(long int));
        if (R[i] == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
    }
    printf("\nn: %.0d\n",n);
    printf("\nTenta estimar a ordem r=ord(x,N) ou multiplo ou divisor dela para extrair os fatores de N\n");
    for (int i=0;i<n;i++){
        for(int j =-1; j<2; j++){
            printf(". ");
            printf("\nresult: %.0f\n",result[i]+j);
            long int *l = FracCont((result[i]+j), q, N, &taml);
            printf("tamanho de l: %d\n", taml);
            printf("[");
            for(int o=0; o<taml;o++){
                printf("%ld ", l[o]);
            }
            printf("]\n");
            long int **t = Frac(l, taml);
            R[i][j+1]=(long int)t[1][0];
            printf("\nt: %ld\n",t[1][0]);
        }
    }
    return R;
}
int main(){
    int x = 2;
    double r =308;
    double q = 1024*1024;
    double N = 23*29;
    int tam=15;
    double result[15] =  {1, 74898, 779624, 953252, 1031555, 289380, 1, 105539, 422155, 255336, 1021340, 674085, 1, 1, 735365 };
    long int **R = EstimaOrdem(r, result, q, N, tam);
   printf("\nR: \n");
    printf("[");
    for(int i=0; i<tam;i++){
        printf("[");
        for(int j=0; j<3; j++){
            printf("%ld ",R[i][j]);
        }
        printf("] ");
        
    }
    printf("]\n");

    return 0;
}