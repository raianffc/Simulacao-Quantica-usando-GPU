#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include<stdbool.h>
#include<math.h>
#include<time.h>
#include<complex.h>
#include <fftw3.h>
#include <unistd.h> 

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

int main(){
    double p1 = 37;
    double p2 = 41;
    double N  = p1 * p2; //N nao precisa ser semi-primo
    double x  = 2;
    double r  = 0;
    double q  = (int)pow(2, 20);//2**20
    int n  = 15; // quantidade de valores medidos 
    float *Soma;
    float *P;
    int tamSoma_P;
    double *Z; 
    float **R;

    return 0;
}