#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include<stdbool.h>
#include<math.h>
#include<time.h>
#include<complex.h>
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

unsigned long* FracCont(unsigned long x, unsigned long q, unsigned long N, int *tamL) {
    int tam = 1;
    double x_inic = (double)x;
    double x_double = (double)x;
    unsigned long *L;
    L = (unsigned long*)malloc(tam * sizeof(unsigned long));
    if (L == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    if (x == 1) {
        tam = 2;
        L = (unsigned long*)realloc(L, tam * sizeof(unsigned long));
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
    x_double = x_double / (double)q;
    int i = 0;
    double max = log(x_inic) / log(1.6);
    do {
        int c = (int) x_double;
        if (tam != 1) {
            L = (unsigned long*)realloc(L, (tam + 1) * sizeof(unsigned long));
            if (L == NULL) {
                printf("Erro na alocacao de memoria.");
                exit(1);
            }
            L[tam] = (unsigned long) c;
        } else {
            L[0] = (unsigned long) c;
        }
        x_double = x_double - (double)c;
        if (x_double >= 0) {
            x_double = 1 / x_double;
            if (x_double > x_inic || x_double > N) {
                x_double = 0;
            }
        } else {
            x_double = 0;
        }
        i++;
        tam++;
    } while (x_double > 0 && i < max);
    *tamL = tam - 1;
    return L;
}

unsigned long** Frac(unsigned long *L, int *tamL) {
    int tamF = *tamL;
    if (*tamL == 0) {
        *tamL = 1;
        L = (unsigned long*)realloc(L, (*tamL) * sizeof(unsigned long));
        if (L == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
        L[0] = 1;
    }
    unsigned long **F = (unsigned long**)malloc(2 * sizeof(unsigned long *));
    if (F == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for (int i = 0; i < 2; i++) {
        F[i] = (unsigned long*)malloc((tamF) * sizeof(unsigned long));
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
        F[0][i+1]=(unsigned long)num;
        F[1][i+1]=(unsigned long)den;
    }
    F[0][0]=(tamF)-1;
    printf("print F[0][0]: %lu", F[0][0]);

    int total = (int)F[1][2];
    for(int j=3;j<tamF;j++){
        total = mmc(total,(int)F[1][j]);
    }
    F[1][0]=(unsigned long)total;
    if(F[1][0]==0 && *tamL==1){
        printf("\nFracao inexistente: %lu \n", L[0]);
    }
    free(L);
    return F;
}

unsigned long **EstimaOrdem(unsigned long r,unsigned long *result,unsigned long q, unsigned long N, int n){
    for(int j =0; j<n; j++){
        printf(" Result[%d]: %lu",j, result[j]);
    }
    int mult = 0;
    unsigned long **R;
    int taml = 1;
    R=(unsigned long**)malloc(n*sizeof(unsigned long*));
    if (R == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for(int i=0;i<n;i++){
        R[i]=(unsigned long*)malloc(3*sizeof(unsigned long));
        if (R[i] == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
    }
    printf("\nTenta estimar a ordem r=ord(x,N) ou multiplo ou divisor dela para extrair os fatores de N\n");
    for (int i=0;i<n;i++){
        for(int j =-1; j<2; j++){
            printf(".");
            unsigned long *l = FracCont(result[i]+j,q,N, &taml);
            unsigned long **t = Frac(l, &taml);
            R[i][j+1]=t[1][0];
            free(t);
        }
        printf("\n");
    }
    return R;
}

int main(){
    unsigned long p1 = 37;
    unsigned long p2 = 41;
    unsigned long N  = p1 * p2; //N nao precisa ser semi-primo
    unsigned long x  = 2;
    unsigned long r  = 0;
    unsigned long q  = (int)pow(2, 20);//2**20
    int n  = 15; // quantidade de valores medidos 
    float *Soma;
    float *P;
    int tamSoma_P;
    unsigned long *Z; 
    unsigned long **R;
    unsigned long result [15]= {343700, 774782, 570891, 367001, 93206, 1042750, 425256, 1, 757304, 1042750, 343700, 1042750, 297096, 75730, 413604};
    R = EstimaOrdem(r, result, q, N, n);
    for(int i=0; i<15; i++){
        for(int j =0; j<3; j++){
            printf(" R[%d][%d]: %lu", i,j, R[i][j]);
        }
        printf("\n");
    }
    return 0;
}