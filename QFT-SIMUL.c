#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
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
            printf("Erro na alocacao de memoria.");
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
    printf("total: %d\n", total);
    F[1][0]=(double)total;
    if(F[1][0]==0 && *tamL==1){
        printf("\nFracao inexistente: %f \n", L[0]);
    }
    free(L);
    return F;
}
/*
float Prepara(float N, float x, float r, float q, float *P, float *Soma){
    r=0;
    q=1;
    int i;
    float s;
 //colocar o valor de r caso seja conhecido (evita o calculo abaixo)
    int tamN = (int)(log2(N));
    float q1 = pow(2,(2*tamN));    // este � o valor ideal segundo Shor. N�o � usado no programa. Serve apenas de refer�ncia
    print("valor ideal para q:",q1);
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
        printf("Ordem r n�o informada. Ordem r calculada: %f",r);
    }
    else{
        printf("Ordem r informada: %f",r);
    }
    printf("Criando Z...");
    float *Z;
    Z= malloc(q*sizeof(float));
    for(int i=0;i<q;i++){
        Z[i]=0;
    }
    int j = 1;
    int cont =0;
    while (j <= q){
        Z[j] = 1;
        j += r;
        cont += 1;  //# em principio, nao e usado. Mas poderia ser usado para normalizar o vetor de estado
    }
//   print(cont)
    //h=cupy.fft.fft(Z)
    //fft_prob(Z,h,Z)
//    mostraQFT(Z)

    printf("Criando Soma com probabilidade acumulada");
    P=malloc((2*r+1)*sizeof(float));
    for(int i=0;i<(2*r+1);i++){
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
        Soma[2*i+1] = total;float Prepara(float N, float x, float r, float q, float *P, float *Soma){
    r=0;
    q=1;
    int i;
    float s;
 //colocar o valor de r caso seja conhecido (evita o calculo abaixo)
    int tamN = (int)(log2(N));
    float q1 = pow(2,(2*tamN));    // este � o valor ideal segundo Shor. N�o � usado no programa. Serve apenas de refer�ncia
    print("valor ideal para q:",q1);
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
        printf("Ordem r n�o informada. Ordem r calculada: %f",r);
    }
    else{
        printf("Ordem r informada: %f",r);
    }
    printf("Criando Z...");
    float *Z;
    Z= malloc(q*sizeof(float));
    for(int i=0;i<q;i++){
        Z[i]=0;
    }
    int j = 1;
    int cont =0;
    while (j <= q){
        Z[j] = 1;
        j += r;
        cont += 1;  //# em principio, nao e usado. Mas poderia ser usado para normalizar o vetor de estado
    }
//   print(cont)
    //h=cupy.fft.fft(Z)
    //fft_prob(Z,h,Z)
//    mostraQFT(Z)

    printf("Criando Soma com probabilidade acumulada");
    P=malloc((2*r+1)*sizeof(float));
    for(int i=0;i<(2*r+1);i++){
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
    P[2*r-1]=-1*(int)(random()*q)
    Soma[((int)2*r-1)]=1;
    printf("Probabilidade Acumulada (dos picos): %f",total);

    return r;
}
    }
    P[2*r-1]=-1*(int)(random()*q)
    Soma[((int)2*r-1)]=1;
    printf("Probabilidade Acumulada (dos picos): %f",total);

    return r;
}*/
double *Simula(float *Soma,float *P, int n, int *tamSoma){
    int i;
    printf("\nSimula medicao:\n");
    double *result;
    result=(double*)malloc(n*sizeof(double));
    for(i = 0; i<n; i++){
        result[i]=0;
    }
     for(i = 0; i<n; i++){
        float m=(rand()/(float)RAND_MAX);
        result[i] = buscabin(Soma,P,m, *tamSoma);
        if (result[i]==0)
            result[i] = 1;
     }
    return result;
}
float **EstimaOrdem(double r,double *result,double q, double N, int *tamresult){
    int mult = 0;
    float **R;
    int taml = 1;
    R=malloc(*tamresult*sizeof(float*));
    for(int i=0;i<*tamresult;i++){
        R[i]=malloc(3*sizeof(float));
    }
    printf("Tenta estimar a ordem r=ord(x,N) ou multiplo ou divisor dela para extrair os fatores de N\n");
    for (int i=0;i<*tamresult;i++){
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
float **EstimaFator(double N, double x,float **R, int tam){
    printf("Procura um multiplo da ordem ou um divisor que distiga um fator nao trivial.\n");
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
            pot=(((int)(pow(x1,R1)))%N);
             if (pot==1){
                sucesso[j]=1;
                printf("multiplo de r\n");
            }
            else{
                int d=mdc(pot-1,N);
                if (d > 1){
                    sucesso[j]=2;
                    printf("fator: %d\n",d);
                }
            }
            potTotal = potTotal*pot;
            if(potTotal%N==1){
                sucesso[j]=sucesso[j]+4;
                printf("fatores do múltiplo da ordem %d %d\n",pot,(potTotal/pot));
            }
            else{
                int d=mdc((potTotal-1),N);
                if (d > 1){
                    sucesso[j]=sucesso[j]+8;
                    printf("valor de sucesso: %d\n",sucesso[j]);
                }
            }
            potTotal = potTotal%N;
            //printf("potTotal: %d", potTotal);
            Sucesso[i][j]=sucesso[j];
        }

    }

    return Sucesso;
}
/*
int* Fatores(N,x,R,S){
    fat=set()
    for i in range(len(S)):
        for j in range(3):
            if S[i][j]==1: # o valor � um m�ltiplo da ordem
                print('testa um divisor da ordem')
                if R[i][j]%2==0:
                    print('Teste de Shor...')
                    print(gcd(pow(x,int(R[i][j]//2),N)-1,N))
                if R[i][j]%3==0:
                    print('Teste com 3 ...')
                    print(gcd(pow(x,int(R[i][j]//3),)-1,N))
            elif S[i][j]==2:   # o valor � um divisor da ordem que encontra um fator
                f=gcd(pow(x,int(R[i][j]),N)-1,N)
                fat.add(f)
                fat.add(N//f)
    return fat
}
*/
int main(){
    int p1 = 5;
    int p2 = 3;
    int N  = p1 * p2; //N n�o precisa ser semi-primo
    int x  = 2;
    int r  = 0;
    int q  = 1024*1024;//2**20
    int n  = 15;// quantidade de valores medidos
    int tamS=1;
    int tamR;

    /*Soma,P,r=Prepara(N,x,r,q);
     #print("sum",Soma)*/
    float *result=Simula(Soma, P, n);
    printf("q= %d\n",q);
    printf("\nResultados medidos na rotina quantica do QOFA:\n");
    printf("[");
    for(int i=0; i<n;i++){
        printf("%.0f ", result[i]);
    }
    printf("]\n");
    int **R = EstimaOrdem(r, result, &tamR);
    /*int **S = EstimaFator(N, x, R, tam);
    printf("[");
   /* for(int i=0; i<n;i++){
        for(int j=0; j<n;j++)
        printf("%d ", S[i][j]);
    }
    printf("]\n");*/

    //int fat = Fatores(N,x,R,S);


    return 0;
}
