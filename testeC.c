#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>
#define PI 3.14159265358979323846

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

int main(){
    double p1 = 29;
    double p2 = 31;
    double N  = p1 * p2; //N n�o precisa ser semi-primo
    double x  = 2;
    double r  = 0;
    double q  = 1024*1024;//2^20
    int n  = 15;// quantidade de valores medidos 
    double *Z;
    float *Soma;
    float *P;
    int tamSoma_P; // variavel P e Soma tem o mesmo tamanho
    
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
    
    return 0;
}
