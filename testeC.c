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
double *Prepara(double N, double x, double *r, double q){
    int tamN = (int)log2(N);
    double q1 = 1 << (2 * tamN);  // este é o valor ideal segundo Shor. Não é usado no programa. Serve apenas de referência
    printf("Valor ideal para q: %.0f\n", q1);

    if(q<N){
        q=1 << (tamN + 4);
    }
    if (*r == 0) {
        int s = x;
        int i = 1;
        while (s > 1) {
            s = (int)(s*x)%((int)N);
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
    fftw_complex *Y = (fftw_complex *)fftw_malloc(q * sizeof(fftw_complex));
    // Otimizar essa parte daqui...
    for (int i = 0; i < q; i++) {
        Y[i] = 0.0 + 0.0 * I;
    }
    int j = 1;
    int cont = 0;
    while (j <= q) {
        Y[j] = 1.0 + 0.0 * I;
        j += *r;
        cont++;
    }
    //Até aqui...
    printf("Calculando FFT...\n");
    fftw_plan plan = fftw_plan_dft_1d(q, Y, Y, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);

    printf("Calculando probabilidades...\n");
    double sum_Z = 0;
    for (int i = 0; i < q; i++) {
        Z[i] = (double)cabs((creal(Y[i])*creal(Y[i]))+(cimag(Y[i])*cimag(Y[i])));
        sum_Z += Z[i];
        //printf("%f --- %f + %fi\n", creal(Z[i]), creal(Y[i]), cimag(Y[i]));
    }
    double temp=0;
    for (int i = 0; i < q; i++) {
        Z[i] = Z[i]/sum_Z;
        temp += Z[i];
    }
    sum_Z = temp;
    printf("Soma das probabilidades: %.20f\nCriando Soma com probabilidade acumulada...\n", sum_Z);
    
    fftw_destroy_plan(plan);
    fftw_free(Y);
    
    return Z;
}
float *Soma_P (double r, double q, float *P, float *Soma, double *Z){
    double k = (q / r);
    printf("Calculando Somas...\n");
    double total = 0;
    for (int i = 0; i <r; i++) {
        double pos = (i * k);
        P[2 * i] = (float)pos;
        total += Z[(int)pos];
        Soma[2*i]= total;
        P[2*i+1]= pos+1;
        total = total + Z[((int)pos)+1];
        Soma[2*i+1] = total;
    }
    P[2*(((int)r)-1)]=-1*(int)((random()%101)*q);
    Soma[2*(((int)r)-1)]=1;

    printf("Probabilidade Acumulada (dos picos): %f\n", total);
    
    return Soma;
}
double *Simula(float *Soma,float *P, int tamResult, int tamSoma){
    srand(time(NULL));
    printf("\nSimula medicao:\n");
    double *result;
    result=(double*)malloc(tamResult*sizeof(double));
    if (result == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for(int i = 0; i<tamResult; i++){
        result[i]=0;
    }
     for(int i = 0; i<tamResult; i++){
        int random_integer = rand();

        // Normaliza o número para estar entre 0 e 1, incluindo valores decimais
        double m = (double)random_integer / (double)RAND_MAX;
        result[i] = (double)buscabin(Soma,P,m, tamSoma);
        if (result[i]==0){
            result[i] = 1;
        }
    }
    return result;
}
float **EstimaOrdem(double r,double *result,double q, double N, int n){
    int mult = 0;
    float **R;
    int taml = 1;
    R=(float**)malloc(n*sizeof(float*));
    if (R == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for(int i=0;i<n;i++){
        R[i]=(float*)malloc(3*sizeof(float));
        if (R[i] == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
    }
    printf("\nTenta estimar a ordem r=ord(x,N) ou multiplo ou divisor dela para extrair os fatores de N\n");
    for (int i=0;i<n;i++){
        for(double j =-1; j<2; j++){
            printf(".");
            double *l = FracCont(result[i]+j,q,N, &taml);
            double **t = Frac(l, &taml);
            int k = (float)j;
            R[i][k+1]=(float)t[1][0];
            free(t);
        }
    }
    return R;
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

    Z = Prepara(N, x, &r, q);
    
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
    
    double *result;
    result = Simula(Soma, P, n, tamSoma_P);
    printf("\nq= %.0f\n",q);
    printf("\nResultados medidos na rotina quantica do QOFA:\n");
    printf("[");
    for(int i=0; i<n;i++){
        printf("%.0f ", result[i]);
    }
    printf("]\n");
    R= EstimaOrdem(r, result, q, N, n);
    free(R);
    return 0;
}