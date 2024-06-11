#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include<stdbool.h>
#include<math.h>
#include<time.h>
#include<complex.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

int mdc(int num1, int num2) {
    int resto;
    do {
        resto = num1 % num2;
        num1 = num2;
        num2 = resto;
    } while (resto != 0);
    print("testando");
    return num1;
}

int mmc(int num1, int num2) {
    int resto, a;
    if(num2==0) return num1;
    a = mdc(num1,num2);
    return (num1 * num2) / a;

}mpz_t *FracCont(int *tamL, double x, int q, int N) {
    int tam = 1;
    double x_inic = x;
    mpz_t *L;

    // Alocação inicial
    L = (mpz_t *)malloc(tam * sizeof(mpz_t));
    if (L == NULL) {
        printf("Erro na alocacao de memoria.\n");
        exit(1);
    }

    mpz_init_set_ui(L[0], 0); // Inicializar a primeira posição de L

    if (x == 1) {
        tam = 2;
        L = (mpz_t *)realloc(L, tam * sizeof(mpz_t));
        if (L == NULL) {
            printf("Erro na alocacao de memoria.\n");
            exit(1);
        }
        mpz_init_set_ui(L[0], 1);
        mpz_init_set_ui(L[1], 0);
        x = 0;
    } else {
        mpz_set_si(L[0], -1);
    }

    x = x / q;
    int i = 0;
    double max = log(x_inic) / log(1.6);

    do {
        int c = (int)x;
        L = (mpz_t *)realloc(L, (tam + 1) * sizeof(mpz_t));
        if (L == NULL) {
            printf("Erro na alocacao de memoria.\n");
            exit(1);
        }
        mpz_init_set_ui(L[tam], c);

        x = x - c;
        if (x >= 0.00001) {  // Use um epsilon pequeno para evitar erros de precisão
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

    *tamL = tam;

    printf("tamanho de L = %d\ntamanho de tam: %d\n", *tamL, tam);
    for (int i = 0; i < tam; i++) {
        gmp_printf("L[%d]: %Zd\n", i, L[i]);
    }

    return L;
}


// Função Frac usando GMP
mpz_t **Frac(mpz_t *L, int *tamL) {
    mpz_t **F;
    int tamF = *tamL;

    if (*tamL == 0) {
        *tamL = 1;
        mpz_realloc2(*L, (*tamL) * sizeof(mpz_t));
        if (L == NULL) {
            printf("Erro na alocacao de memoria.\n");
            exit(1);
        }
        mpz_init_set_ui(L[0], 1);
    }

    F = (mpz_t **)malloc(2 * sizeof(mpz_t *));
    for (int i = 0; i < 2; i++) {
        F[i] = (mpz_t *)malloc(tamF * sizeof(mpz_t));
        for (int j = 0; j < tamF; j++) {
            mpz_init(F[i][j]);
        }
    }

    mpz_t num, den, temp;
    mpz_init(num);
    mpz_init(den);
    mpz_init(temp);

    for (int i = 0; i < tamF; i++) {
        mpz_set_ui(num, 1);
        mpz_set(den, L[tamF - i - 1]);
        for (int j = tamF - i - 2; j >= 0; j--) {
            mpz_set(temp, num);
            mpz_set(num, den);
            mpz_mul(den, L[j], den);
            mpz_add(den, den, temp);
        }
        mpz_set(F[0][i], num);
        mpz_set(F[1][i], den);
    }

    mpz_t total;
    mpz_init(total);
    mpz_set(total, F[1][0]);

    for (int j = 1; j < tamF; j++) {
        mpz_lcm(total, total, F[1][j]);
    }

    mpz_set(F[1][0], total);

    mpz_clear(num);
    mpz_clear(den);
    mpz_clear(temp);
    mpz_clear(total);

    return F;
}

mpz_t **EstimaOrdem(unsigned long r, double *result, unsigned long q, unsigned long N, int n) {
    mpz_t **R;

    // Alocar memória para R usando GMP
    R = malloc(n * sizeof(mpz_t *));
    if (R == NULL) {
        printf("Erro na alocação de memória.\n");
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        R[i] = malloc(3 * sizeof(mpz_t));
        if (R[i] == NULL) {
            printf("Erro na alocação de memória.\n");
            exit(1);
        }
        for (int j = 0; j < 3; j++) {
            mpz_init(R[i][j]);
        }
    }

    printf("\nTentando estimar a ordem r=ord(x,N) ou múltiplo ou divisor dela para extrair os fatores de N\n");

    for (int i = 0; i < n; i++) {
        for (int j = -1; j <= 1; j++) {
            printf(".");

            // Chamar FracCont com os tipos originais (unsigned long)
            int taml = 0;
            double adjusted_result = result[i] + j; // Ajustando resultado
            mpz_t *frac_cont_result = FracCont(&taml, adjusted_result, q, N);

            // Chamar Frac com os tipos originais (unsigned long)
            mpz_t **t = Frac(frac_cont_result, &taml);

            gmp_printf("t[1][0]: %Zd\n", t[1][0]);

            // Converter o resultado de Frac para mpz_t para armazenar em R
            mpz_set(R[i][j + 1], t[1][0]);

            // Liberação da memória para frac_cont_result
            for (int k = 0; k < taml; k++) {
                mpz_clear(frac_cont_result[k]);
            }
            free(frac_cont_result);

            // Liberação da memória para t
            for (int k = 0; k < 2; k++) {
                for (int l = 0; l < taml; l++) {
                    mpz_clear(t[k][l]);
                }
                free(t[k]);
            }
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
    mpz_t **R;
    double result [15]= {343700, 774782, 570891, 367001, 93206, 1042750, 425256, 1, 757304, 1042750, 343700, 1042750, 297096, 75730, 413604};
    R = EstimaOrdem(r, result, q, N, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++) {
            gmp_printf("R[%d][%d]: %Zd\n", i, j, R[i][j]);
            mpz_clear(R[i][j]);
        }
        free(R[i]);
    }
    free(R);
    return 0;
}
