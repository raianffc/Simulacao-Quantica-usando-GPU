#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>
#define PI 3.14159265358979323846

int* Fatores(double N, double x, double **R, int **S, int num_s, int num_j, int *num_fatores) {
    int *fat = (int *)malloc((num_s * num_j * 2) * sizeof(int));
    int count = 0;

    for (int i = 0; i < num_s; i++) {
        for (int j = 0; j < num_j; j++) {
            if (S[i][j] == 1) {
                if ((int)R[i][j] % 2 == 0) {
                    int f = gcd((int)(pow(x, (int)(R[i][j] / 2)) - 1), (int)N);
                    fat[count++] = f;
                }
                if ((int)R[i][j] % 3 == 0) {
                    int f = gcd((int)(pow(x, (int)(R[i][j] / 3)) - 1), (int)N);
                    fat[count++] = f;
                }
            } else if (S[i][j] == 2) {
                int f = gcd((int)pow(x, (int)R[i][j]) - 1, (int)N);
                fat[count++] = f;
                fat[count++] = (int)N / f;
            }
        }
    }

    *num_fatores = count;
    return fat;
}
int* Fatores( double N, double x, float **R, float **S, int *tamR, int *tamS, int *k){
    int *fat;
    fat = (float*)malloc((*k)*sizeof(int));
    int f;
    if(fat==NULL){ 
    	printf("\nerror\n");
    	exit(1);
    }
    for (int i=0; i<*tamS; i++){
        for(int j=0; j<3; j++){
            if (S[i][j]==1){ // o valor � um multiplo da ordem
                printf("\ntesta um divisor da ordem\n");
                if(((int)R[i][j])%2==0){
                    printf("\nTeste de Shor...\n");
                    printf(mdc(pow(x,((int)R[i][j]/2),((int)N))-1,N)); //aqui
                }
                if(((int)R[i][j])%3==0){
                    print("Teste com 3 ...");
                    print(mdc(pow(x,((int)R[i][j]/2),)-1,N));//aqui
                }
            }
            else if (S[i][j]==2){   // o valor e um divisor da ordem que encontra um fator
                f=mdc(pow(x,((int)R[i][j]/2),N)-1,N);//aqui
                //Minha outra dúvida seria como fazer essa lista fat ser imutavel que nem o Python faz usando set que é uma coleção
		        fat[(*k)-1]=f;
		        *k=(*k)+2;
		        fat=realloc(fat, (*k) * sizeof(int));
                fat[(*k)-2]=(int)(N/f);
     	 	}
        }
    }
    return fat;
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
