#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include<stdbool.h>
#include<math.h>
#include<time.h>
#include<complex.h>
#include <unistd.h> 
#include <gmp.h>
#include <stdlib.h>
#include <assert.h>

int mdc(unsigned int num1, unsigned int num2) {
    unsigned int resto;
    do {
        resto = num1 % num2;
        num1 = num2;
        num2 = resto;
    } while (resto != 0);
    return num1;
}

int mmc(unsigned int num1, unsigned int num2) {
    unsigned int resto, a;
    if(num2==0) return num1;
    a = mdc(num1,num2);
    return (num1 * num2) / a;
    
}

int* FracCont(double x, int q, int N, int *tamL) {
    int tam = 1;
    double x_inic = x;
    int *L;
    L = (int*)malloc(tam * sizeof(int));
    if (L == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    if (x == 1) {
        tam = 2;
        L = (int*)realloc(L, tam * sizeof(int));
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
    printf("\n\n\nmax: %f\n\n\n", max);
    do {
        int c = (int) x;
        if (tam != 1) {
            L = (int*)realloc(L, (tam + 1) * sizeof(int));
            if (L == NULL) {
                printf("Erro na alocacao de memoria.");
                exit(1);
            }
            L[tam] = (int) c;
        } else {
            L[0] = (int) c;
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
    for(int i =0; i<tam; i++){
        printf("L[%d]: %d\n", i,L[i]);
    }
    return L;
}

int** Frac(int *L, int *tamL) {
    int tamF = *tamL;
    if (*tamL == 0) {
        *tamL = 1;
        L = (int*)realloc(L, (*tamL) * sizeof(int));
        if (L == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
        L[0] = 1;
    }
    int **F = (int**)malloc(2 * sizeof(int *));
    if (F == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for (int i = 0; i < 2; i++) {
        F[i] = (int*)malloc((tamF) * sizeof(int));
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
        F[0][i+1]=(int)num;
        F[1][i+1]=(int)den;
    }
    F[0][0]=(tamF)-1;

    int total = (int)F[1][2];
    for(int j=3;j<tamF;j++){
        total = mmc(total,(int)F[1][j]);
    }
    
    F[1][0]=(int)total;
    if(F[1][0]==0 && *tamL==1){
        printf("\nFracao inexistente: %d \n", L[0]);
    }
    free(L);
    return F;
}
int **EstimaOrdem(int r,double *result,int q, int N, int n){
    int mult = 0;
    int **R;
    int taml = 1;
    R=(int**)malloc(n*sizeof(int*));
    if (R == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for(int i=0;i<n;i++){
        R[i]=(int*)malloc(3*sizeof(int));
        if (R[i] == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
    }
    printf("\nTenta estimar a ordem r=ord(x,N) ou multiplo ou divisor dela para extrair os fatores de N\n");
    for (int i=0;i<n;i++){
        for(int j =-1; j<2; j++){
            printf(". ");
            int *l = FracCont(result[i]+j,q,N, &taml);
            int **t = Frac(l, &taml);
            printf(" t[1][0]: %d", t[1][0]);
            printf("\n");
            R[i][j+1]=t[1][0];
            free(t);
        }
    }
    return R;
}

int main(){

    int p1 = 37;
    int p2 = 41;
    int N  = p1 * p2; //N nao precisa ser semi-primo
    double x  = 2; //Tem que ser do tipo double
    int r  = 0;
    int q  = (int)pow(2, 20);//2**20
    int n  = 15; // quantidade de valores medidos 
    float *Soma;
    float *P;
    int tamSoma_P;
    int **R;
    float **S;
    int *fat;
    int tamFat;
    double *Z; 
    //Result tem que ser double
    double result [15]= {343700, 774782, 570891, 367001, 93206, 1042750, 425256, 1, 757304, 1042750, 343700, 1042750, 297096, 75730, 413604};
    R = EstimaOrdem(r, result, q, N, n);

    /*Teste do gmp*/
     mpz_t matrix[15][3];
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 3; j++) {
            mpz_init(matrix[i][j]);
        }
    }
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 3; j++) {
            mpz_set_ui(matrix[i][j], R[i][j]); // Atribuir valores de exemplo
        }
    }
     for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 3; j++) {
            gmp_printf("matrix[%d][%d] = %Zd\n", i, j, matrix[i][j]);
        }
    }
    // Limpar a memória usada pelos elementos da matriz
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 3; j++) {
            mpz_clear(matrix[i][j]);
        }
    }
    return 0;
}
