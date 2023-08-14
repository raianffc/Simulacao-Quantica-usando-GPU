#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>
#define PI 3.14159265358979323846

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


int main(){
    double p1 = 29;
    double p2 = 31;
    double N  = p1 * p2; //N n�o precisa ser semi-primo
    double x  = 2;
    double r  = 0;
    double q  = 1024*1024;//2^20
    int n  = 15;// quantidade de valores medidos 
    float *Soma;
    int tamSoma;
    float *P; 
    int tamP;
    r = Prepara(N, x, r, q, P, Soma, &tamP, &tamSoma);
    printf("%f", r);
    return 0;
}
