/*
 * main.c
 *
 *  Created on: 1 de mai de 2017
 *  	Author: cx3d
 */

#include <math.h>
#include <stdio.h>

#define FALSE 0
#define TRUE 1
#define MAXITER 100

typedef double (*funcao_t)(double);

double g1(double x) {
    return (x * x + 2.0) / 3.0;
}

double g2(double x) {
    return sqrt(3.0 * x - 2.0);
}

double g3(double x) {
    return 3.0 - 2.0 / x;
}

double g4(double x) {
    return (x * x - 2.0) / (2.0 * x - 3.0);
}

funcao_t funcoes[4] = { g1, g2, g3, g4 };

int iteracao_ponto_fixo(double x0, int idx, double tol, double* x) {
    funcao_t g = funcoes[idx];
    *x = g(x0);
    return fabs((*x) - x0) < tol;
}

int main(int argc, char** argv) {
    double tol = 1.0e-6;

    double xstart;
    printf("Digite estimativa inicial de x:");
    scanf("%lf", &xstart);

    for (int idx = 0; idx < 4; idx++) {
   	 double x;
   	 double x0 = xstart;
   	 int convergiu = FALSE;
   	 int i;
   	 for (i = 0; (!convergiu) & (i < MAXITER); i++) {
   		 convergiu = iteracao_ponto_fixo(x0, idx, tol, &x);
   		 x0 = x;
   	 }

   	 if (i == MAXITER)
   		 printf("Não convergiu para a função de ponto fixo %d.\n", (idx + 1));
   	 else
   		 printf("Solução encontrada para a função %d: %f em %d iterações.\n", (idx + 1), x, i);
    }

    return 0;

}


