#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
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
double Prepara(double N, double x, double r, double q, float *Soma,float *P, int *tamS, int *tamP){
    r=0;
    q=1;
    int i;
    float s;
 //colocar o valor de r caso seja conhecido (evita o calculo abaixo)
    int tamN = (int)(log2(N));
    float q1 = pow(2,(2*tamN));    // este � o valor ideal segundo Shor. N�o � usado no programa. Serve apenas de refer�ncia
    printf("\nvalor ideal para q: %f\n",q1);
    if(q < N)
        q = 1 << (tamN+4); // bitwise deslocamento s esquerda
    if (r==0){
        s=x;
        i=1;

        while (s > 1){
            s = ((int)(s*x))%((int)N);
            i++;
        }
        r = i;
        printf("\nOrdem r nao informada. Ordem r calculada: %f\n",r);
    }
    else{
        printf("\nOrdem r informada: %f\n",r);
    }
    printf("\nCriando Z...\n");
    double complex *Z;
    Z= malloc(q*sizeof(double complex));
    for(int i=0;i<q;i++){
        Z[i]=0;
    }
    int j = 1;
    int cont =0;
    while (j <= q){
        Z[j] = 1;
        j = j +((int)r);
        cont ++;  //# em principio, nao e usado. Mas poderia ser usado para normalizar o vetor de estado
    }
    /*printf("[");
    for(int i=0; i<q;i++){
        printf("%.1f \n", creal(Z[i]));
    }
    printf("]\n");*/
//   print(cont)
    printf("Calculando FFT...");
    double soma=0;
    ifft(Z, (int)N);
    /*printf("[");
    for(int i=0; i<q;i++){
        printf("%.1f \n", creal(Z[i]));
    }
    printf("]\n");*/

    printf("Calculando probabilidades...");
    for(int i=0; i<q; i++){
        Z[i] = abs(Z[i]*Z[i]);
        soma = soma + Z[i];
    }
    printf("\nsoma = %f\n", soma);
    for(int i=0; i<q; i++){
        Z[i] = Z[i]/soma; //normalizando vetor ??
    }
    /*printf("[");
    for(int i=0; i<q;i++){
        printf("%.1f \n", creal(Z[i]));
    }
    printf("]\n");*/

//    mostraQFT(Z)

    printf("Criando Soma com probabilidade acumulada");
    P=malloc((2*(r+1))*sizeof(float));
    for(int i=0;i<(2*(r+1));i++){
        P[i]=0;
    }
    Soma=malloc((2*r+1)*sizeof(float));
    for(int i=0;i<(2*r+1);i++){
        Soma[i]=0;
    }
    float k = q/r;
    printf("Calculando Somas...");
    float total = 0;
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
/*
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
*/
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
    int k=1;
    
    
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

    //fat = Fatores(N, x, R, S, n, n, &k);
    //print("Tamanho de k e %d", k);
    
    return 0;
}
