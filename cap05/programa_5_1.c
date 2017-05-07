/*
  * main.c
  *
  *  Created on: 1 de mai de 2017
  *          Author: cx3d
  */

#include <math.h>
#include <stdio.h>

#define FALSE 0
#define TRUE 1
#define MAXITER 100

typedef double (*funcao_t)(double);

double g1(double x) {
	return x * x - 2.0;
}

double g2(double x) {
	return sqrt(x + 2.0);
}

double g3(double x) {
	return 1.0 + 2.0 / x;
}

double g4(double x) {
	return (x * x + 2.0) / (2.0 * x - 1.0);
}

double C[4] = { 4.0, 0.25, 0.5, 2.0 / 3.0 };

funcao_t funcoes[4] = { g1, g2, g3, g4 };

int iteracao_ponto_fixo(double x0, funcao_t g, double tol, double* x, double* err) {
	*x = g(x0);
	*err = fabs((*x) - x0);
	return *err < tol;
}

int main(int argc, char** argv) {
	double tol = 1.0e-6;

	double xstart;
	printf("Digite estimativa inicial de x:");
	scanf("%lf", &xstart);

	double r;
	double err, err_prev;

	for (int idx = 0; idx < 4; idx++) {
		printf("\nFuncao %d:\n", (idx + 1));
		printf("Iteracao	Taxa de convergencia:\n");
		double x;
		double x0 = xstart;
		int convergiu = FALSE;
		int i;
		for (i = 0; (!convergiu) & (i < MAXITER); i++) {
			convergiu = iteracao_ponto_fixo(x0, funcoes[idx], tol, &x, &err);
			x0 = x;

			if (i > 0) {
				r = (log(err) - log(C[idx])) / log(err_prev);
				printf("%d\t\t%f\n", i, r);
			}

			err_prev = err;
		}

		if (i == MAXITER)
			printf("Não convergiu para a função de ponto fixo %d.\n", (idx + 1));
		else
			printf("Solução encontrada para a função %d: %f em %d iterações.\n", (idx + 1), x, i);
	}

	return 0;

}
