#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include<stdbool.h>
#include<math.h>
#include<time.h>
#include<complex.h>
#include <unistd.h>
#define PI 3.14159265358979323846

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
    long int a;
    if(num2==0) return num1;
    a = mdc(num1,num2);
    return (num1 * num2) / a;

}

long int* FracCont(double x, double q, double N, int *tamL) {
    int tam = 1;
    if(x<0) x=-1*x;
    double x_inic = x;
    long int *L;
    L = (long int*)malloc(tam * sizeof(long int));
    if (L == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    if (x == 1) {
        tam = 2;
        L = (long int*)realloc(L, tam * sizeof(long int));
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
    double max = (double)log(x_inic) / (double)log(1.6);
    printf("max: %lf\n", max);
    do {
        long int c = (long int)x;
        if (tam != 1) {
            L = (long int*)realloc(L, (tam + 1) * sizeof(long int));
            if (L == NULL) {
                printf("Erro na alocacao de memoria.");
                exit(1);
            }
            L[tam] = c;
        } else {
            L[0] =  c;
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
    printf("max: %d\n", i);
    *tamL = tam - 1;
    return L;
}
/*
long int** Frac(long int*L, int tamL) {
    int tamF = tamL;
    if (tamF == 0) {
        tamF = 1;
        L = (long int*)realloc(L, (tamF) * sizeof(long int));
        if (L == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
        L[0] = 1;
    }
   long int **F = (long int**)malloc(2 * sizeof(long int *));
    if (F == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for (int i = 0; i < 2; i++) {
        F[i] = (long int*)malloc((tamF) * sizeof(long int));
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
    long int num = 1;
    long int den = 1;
    for ( int i=0; i<((tamF)-1);i++){
        num = 1;
        den = L[(tamF)-i-1];
        for (int j=((tamF)-i-2); j>0;j--){
            long int temp=num;
            num = den;
            den = L[j]*den+temp;
        }
        F[0][i+1]= num;
        F[1][i+1]= den;
    }
    F[0][0]=(tamF)-1;

    long int total = F[1][2];
    for(int j=3;j<tamF;j++){
        total = mmc(total, F[1][j]);
    }
    F[1][0]=total;
    if(F[1][0]==0 && tamF==1){
        printf("\nFracao inexistente: %.0f \n", L[0]);
    }
    for (int i=0; i<2; i++){
        for( int j=0; j<tamF; j++){
            printf("F[%d][%d] = %ld\n", i, j, F[i][j]);
        }
    }
    printf("\n\ntotal: %ld \n\n", total);
    free(L);
    return F;
}
long int **EstimaOrdem(double r,long int *result,double q, double N, int n){
    int mult = 0;
    long int **R;
    int taml = 1;
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
    printf("\nTenta estimar a ordem r=ord(x,N) ou multiplo ou divisor dela para extrair os fatores de N\n");
    for (int i=0;i<3;i++){
        for(int j =-1; j<n; j++){
            printf(". ");
            long int *l = FracCont(result[i]+j,q,N, &taml);
            long int **t = Frac(l, taml);
            R[i][j+1]=t[1][0];
            free(t);
        }
    }
    return R;
}*/
int main(){
    double x = 2;
    double q = 1024*1024;
    double N = 47*43;
    double r = 322;
    int n =15;
    int tam;
    long int resultado [15] = {26051, 179105, 696880, 898779, 113976, 709905, 1045319, 719675, 436363, 169335, 280054, 856446, 179105, 807599, 595929};
    long int *result;
    result = (long int*)malloc(n*sizeof(long int));
    for (int i=0; i<n; i++){
        result[i] = resultado[i];
    }
   long int *L = FracCont(x, q, N, &tam);
    for (int i=0; i<n; i++){
        printf("L[%d] = %ld ", i, L[i]);
    }

    /*long int **R = EstimaOrdem(r, result, q, N, n);
    for (int i=0; i<n; i++){
        for( int j=0; j<3; j++){
            printf("R[%d][%d] = %d\n", i, j, R[i][j]);
        }
    }*/

    return 0;
}
