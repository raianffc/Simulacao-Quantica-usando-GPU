#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include<complex.h>
#define PI 3.14159265358979323846

void ifft(double complex *x, int N) {
    if (N <= 1)
        return;

    double complex *even = malloc(N / 2 * sizeof(double complex));
    double complex *odd = malloc(N / 2 * sizeof(double complex));

    for (int i = 0; i < N / 2; i++) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    ifft(even, N / 2);
    ifft(odd, N / 2);

    for (int k = 0; k < N / 2; k++) {
        double complex t = cexp(-I * 2 * PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }

    free(even);
    free(odd);
}

double Prepara(double N, double x, double r, double q, float *Soma,float *P, int *tamS, int *tamP){
    int i;
    double s;
 //colocar o valor de r caso seja conhecido (evita o calculo abaixo)
    int tamN = (int)(log2(N));
    float q1 = pow(2,(2*tamN));    // este eh o valor ideal segundo Shor. Nao eh usado no programa. Serve apenas de refer�ncia
    printf("\nvalor ideal para q: %.0f\n",q1);
    if(q < N){
        q = 1 << (tamN+4); // bitwise deslocamento s esquerda
    }
    if (r==0){
        s=x;
        i=1;
        while (s > 1){
            s = (int)(s*x)%((int)N);
            i++;
        }
        r = i;
        printf("\nOrdem r nao informada. Ordem r calculada: %f\n",r);
    }
    else{
        printf("\nOrdem r informada: %f\n",r);
    }
    printf("\nCriando Z...\n");
    double complex *Y;
    Y= malloc(q*sizeof(double complex));
    for(int i=0;i<q;i++){
        Y[i]=0;
    }
    double *Z;
    Z= malloc(q*sizeof(double));
    int j = 1;
    int cont =0;
    while (j <= q){
        Y[j] = 1;
        j = j +((int)r);
        cont ++;  //# em principio, nao e usado. Mas poderia ser usado para normalizar o vetor de estado
    }
//   print(cont)
    printf("\nCalculando FFT...\n");
    double soma=0;
    ifft(Y, (int)q);
 
    printf("\nCalculando probabilidades...\n");
    for(int i=0; i<q; i++){
        Z[i] = (double)cabsf(csqrt(Y[i]));
        printf("\n%f\n", cimag(Y[i]));
        
        
    }
    soma= soma+Z[0];
    printf("\nSOma = %f\n", soma);
    for(int i=0; i<q; i++){
        soma= soma+Z[i];
    }
    printf("\nsoma = %f\n", soma);
    /*printf("[");
    for(int i=0; i<q;i++){
        printf("%.1f \n", creal(Z[i]));
    }
    printf("]\n");*/


//    mostraQFT(Z)

    printf("\nCriando Soma com probabilidade acumulada\n");
    P=malloc((2*(r+1))*sizeof(float));
    for(int i=0;i<(2*(r+1));i++){
        P[i]=0;
    }
    Soma=malloc((2*r+1)*sizeof(float));
    for(int i=0;i<(2*r+1);i++){
        Soma[i]=0;
    }
    float k = q/r;
    printf("\nCalculando Somas...\n");
    double total = 0;
    for (int i=0;i<r;i++){
        int pos = (int)(i*k);
        P[2*i] = pos;
        total =total+ Z[pos];
        Soma[2*i]= total;
        P[2*i+1]= pos+1;
        total = total + Z[pos+1];
        Soma[2*i+1] = total;
    }
    P[2*(((int)r)-1)]=-1*(int)((random()%101)*q);
    Soma[2*(((int)r)-1)]=1;
    printf("\nProbabilidade Acumulada (dos picos): %f \n",total);

    return r;
    
}
int main(){
    double p1 = 29;
    double p2 = 31;
    double N  = p1 * p2; //N n�o precisa ser semi-primo
    double x  = 1024;
    double r  = 0;
    double q  = 1024*1024;//2**20
    int n  = 15;// quantidade de valores medidos 
    float *Soma;
    int tamS;
    float *P; 
    int tamP;
    r = Prepara(N, x, r, q,Soma, P, &tamS, &tamP);
    printf("%f", r);
    return 0;
}
