#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<complex.h>
#include <fftw3.h>
#define PI 3.14159265358979323846

int mdc(int num1, int num2) {
    int resto;
    do {
        resto = num1 % num2;
        num1 = num2;
        num2 = resto;
    } while (resto != 0);
    return num1;
}

int mmc(int num1, int num2) {
    int resto, a;
    if(num2==0) return num1;
    a = mdc(num1,num2);
    return (num1 * num2) / a;
    
}
float buscabin(float *Soma, float *P, int m, int tam) {
    if( tam == 0)
        return 0;
    else if (tam == 1)
        return P[0];
    else if (tam == 2){
        if (m <= Soma[0]){
            return P[0];
        }
        else{
            return P[1];
        }
    }
    else{
        int meio = tam/2;

        if (m == Soma[meio])
            return P[meio];

        else if (m < Soma[meio])
            return buscabin(Soma,P,m, meio);

        else{
            float r = buscabin(&Soma[meio+1],&P[meio+1],m, (tam-1-meio));
            if(r==-1)
                return -1;

        }
    }
}

double* FracCont(double x, double q, double N, int *tamL) {
    int tam = 1;
    double x_inic = x;
    double *L = malloc(tam * sizeof(double));
    if (L == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    if (x == 1) {
        tam = 2;
        L = realloc(L, tam * sizeof(double));
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
            L = realloc(L, (tam + 1) * sizeof(double));
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
        L = realloc(L, (*tamL) * sizeof(double));
        if (L == NULL) {
            printf("Erro na alocacao de memoria.");
            exit(1);
        }
        L[0] = 1;
    }
    double **F = malloc(2 * sizeof(double *));
    if (F == NULL) {
        printf("Erro na alocacao de memoria.");
        exit(1);
    }
    for (int i = 0; i < 2; i++) {
        F[i] = malloc((tamF) * sizeof(double));
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
        printf("\nFracao inexistente: %f \n", L[0]);
    }
    free(L);
    return F;
}
double Prepara(double N, double x, double r, double q, float *P, float *Soma, int *tamP, int *tamSoma){
    int tamN = (int)log2(N);
    int q1 = 1 << (2 * tamN);  // este é o valor ideal segundo Shor. Não é usado no programa. Serve apenas de referência
    printf("Valor ideal para q: %d\n", q1);

    if(q<N){
        q=1 << (tamN + 4);
    }
    if (r == 0) {
        int s = x;
        int i = 1;
        while (s > 1) {
            s = (int)(s*x)%((int)N);
            i++;
        }
        r = i;
        printf("Ordem r não informada. Ordem r calculada: %.0f\n", r);
    } else {
        printf("Ordem r informada: %f\n", r);
    }

    printf("Criando Z...\n");
    double *Z = (double *)malloc(q * sizeof(double));
    fftw_complex *Y = (fftw_complex *)fftw_malloc(q * sizeof(fftw_complex));
    // Otimizar essa parte daqui...
    for (int i = 0; i < q; i++) {
        Y[i] = 0.0 + 0.0 * I;
    }
    int j = 1;
    int cont = 0;
    while (j <= q) {
        Y[j] = 1.0 + 0.0 * I;
        j += r;
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
    P=malloc((2*(r+1))*sizeof(float));
    /*for(int i=0;i<(2*(r+1));i++){
        P[i]=0;
    }*/
    Soma=malloc((2*r+1)*sizeof(float));
    /*for(int i=0;i<(2*r+1);i++){
        Soma[i]=0;
    }*/
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

    free(Z);
    fftw_destroy_plan(plan);
    fftw_free(Y);
    
    return r;
}
double *Simula(float *Soma,float *P, int n, int *tamSoma){
    int i;
    printf("\nSimula medicao:\n");
    double *result;
    result=(double*)malloc(n*sizeof(double));
    for(i = 0; i<n; i++){
        result[i]=0;
    }
     for(i = 0; i<n; i++){
        float m=((random()%101)/(float)RAND_MAX);
        result[i] = buscabin(Soma,P,m, *tamSoma);
        if (result[i]==0)
            result[i] = 1;
     }
    return result;
}
float **EstimaOrdem(double r,double *result,double q, double N, int n){
    int mult = 0;
    float **R;
    int taml = 1;
    R=malloc(n*sizeof(float*));
    for(int i=0;i<n;i++){
        R[i]=malloc(3*sizeof(float));
    }
    printf("\nTenta estimar a ordem r=ord(x,N) ou multiplo ou divisor dela para extrair os fatores de N\n");
    for (int i=0;i<n;i++){
        for(double j =-1; j<2; j++){
            printf(". ");
            double *l = FracCont(result[i]+j,q,N, &taml);
            double **t = Frac(l, &taml);
            int k = (float)j;
            R[i][k+1]=(float)t[1][0];
            free(t);
        }
    }
    return R;
}

//0 - não é múltiplo e nem divisor da ordem que distinga fator
// 1 - múltiplo da ordem
// 2 - distingue fator
// 4 - mmc múltiplo da ordem (acumulada)
// 8 - mmc distingue fator

float **EstimaFator(double N, double x,float **R, int tam){
    printf("\nProcura um multiplo da ordem ou um divisor que distiga um fator nao trivial.\n");
    float **Sucesso;
    int potTotal=1;
    int pot;
    Sucesso=malloc(tam*sizeof(int));
    for(int i = 0; i<tam; i++){
        Sucesso[i]=malloc(3*sizeof(int));
    }
    for(int i = 0; i<tam; i++){
        float sucesso[3]={0,0,0};
        for(int j=0; j<3; j++){
            double x1=(double)x;
            double R1=(double)R[i][j];
            pot=(((int)(pow(x1,R1)))%((int)N));
             if (pot==1){
                sucesso[j]=1;
                printf("\nmultiplo de r\n");
            }
            else{
                int d=mdc(pot-1,N);
                if (d > 1){
                    sucesso[j]=2;
                    printf("\nfator: %d\n",d);
                }
            }
            potTotal = potTotal*pot;
            if(potTotal%((int)N)==1){
                sucesso[j]=sucesso[j]+4;
                printf("\nfatores do múltiplo da ordem %d %d\n",pot,(potTotal/pot));
            }
            else{
                int d=mdc((potTotal-1),N);
                if (d > 1){
                    sucesso[j]=sucesso[j]+8;
                    printf("\nvalor de sucesso: %f\n",sucesso[j]);
                }
            }
            potTotal = potTotal%((int)N);
            //printf("potTotal: %d", potTotal);
            Sucesso[i][j]=sucesso[j];
        }

    }

    return Sucesso;
}

int* Fatores(double N, double x, float **R, float **S, int num_s, int *num_fatores) {
    int *fat = (int *)malloc((num_s * num_s * 2) * sizeof(int));
    int count = 0;

    for (int i = 0; i < num_s; i++) {
        for (int j = 0; j < 3; j++) {
            if (S[i][j] == 1) {
                if ((int)R[i][j] % 2 == 0) {
                    int f = mdc((int)(pow(x, (int)(R[i][j] / 2)) - 1), (int)N);
                    fat[count++] = f;
                }
                if ((int)R[i][j] % 3 == 0) {
                    int f = mdc((int)(pow(x, (int)(R[i][j] / 3)) - 1), (int)N);
                    fat[count++] = f;
                }
            } else if (S[i][j] == 2) {
                int f = mdc((int)pow(x, (int)R[i][j]) - 1, (int)N);
                fat[count++] = f;
                fat[count++] = (int)N / f;
            }
        }
    }

    *num_fatores = count;
    return fat;
}

int main(){
    double p1 = 29;
    double p2 = 31;
    double N  = p1 * p2; //N n�o precisa ser semi-primo
    double x  = 6440;
    double r  = 0;
    double q  = 1024*1024;//2**20
    int n  = 15;// quantidade de valores medidos 
    float *Soma;
    int tamSoma;
    float *P; 
    int tamP;
    double *result;
    float **R;
    float **S;
    int *fat;
    int Tamfat;
    
    
    r = Prepara(N, x, r, q,Soma, P, &tamSoma, &tamP);
    result = Simula(Soma, P, n, &tamSoma);
    printf("\nq= %f\n",q);
    printf("\nResultados medidos na rotina quantica do QOFA:\n");
    printf("[");
    for(int i=0; i<n;i++){
        printf("%.0f ", result[i]);
    }
    printf("]\n");

    R= EstimaOrdem(r, result, q, N, n);

    S = EstimaFator(N, x, R, n);

    fat = Fatores(N, x, R, S, n, &Tamfat);
    printf("Fatores: \n");
    for(int i=0; i<Tamfat;i++){
        printf("\n%d\n", fat[i]);
    }
    return 0;
}
