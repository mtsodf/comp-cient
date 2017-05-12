#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cblas.h"

#define FALSE 0
#define TRUE 1
#define MAXITER 100

typedef int bool;

/**
 * Imprime um relatorio sucinto sobre a solucao.
 */
void print_solucao(bool convergiu, int numiter, int N, double* x) {
	if (convergiu) {
		printf("Numero de iteracoes: %d\n", numiter);
		printf("Solucao: [%f", x[0]);

		for (int j = 1; j < N; j++) {
			printf(", %f", x[j]);
		}

		printf("]\n");
	} else
		printf("Nao convergiu em %d iteracoes.\n", MAXITER);
}

/*
 * Resolve o sistema linear usando steepest descent.
 * N <- ordem do sistema
 * A <- matriz de coeficientes NxN, representada como um vetor organzado como row-major
 * x <- vetor de estimativa inicial, retornara a solução se o problema convergir
 * b <- segundo membro do sistema linear, retornara o vetor de residuos se o problema convergir
 * tol <- tolerancia para a norma do residuo
 * numiter <- retorna o numero de iteracoes ate a solucao ser atingida
 * Obs.: algoritmo pag 142 do livro do Yousef Saad
 * Link: http://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf
 */
bool solve_steepest_descent(int N, double* A, double* x, double* b, double tol, int* numiter) {
	double* p = (double*) malloc(N * sizeof(double));
	//r <- b - A*x
	cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, -1.0, A, N, x, 1, 1.0, b, 1);

	//p <- A*r
	cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, A, N, b, 1, 0.0, p, 1);

	int convergiu = FALSE;
	int i;
	for (i = 0; (!convergiu) && (i < MAXITER); i++) {
		//alpha = <r, r>/<p, r>
		double alpha1 = cblas_ddot(N, b, 1, b, 1);
		double alpha2 = cblas_ddot(N, p, 1, b, 1);
		double alpha = alpha1 / alpha2;

		//x <- x + alpha*r
		cblas_daxpy(N, alpha, b, 1, x, 1);

		//r <- r - alpha*p
		cblas_daxpy(N, -alpha, p, 1, b, 1);

		double res = cblas_dnrm2(N, b, 1);
		if (res < tol)
			convergiu = TRUE;
		else
			//p <- A*r
			cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, A, N, b, 1, 0.0, p, 1);
	}

	free(p);

	*numiter = i;

	return convergiu;
}

/*
 * Resolve o sistema linear usando gradiente conjugado.
 * N <- ordem do sistema
 * A <- matriz de coeficientes NxN, representada como um vetor organzado como row-major
 * x <- vetor de estimativa inicial, retornara a solução se o problema convergir
 * b <- segundo membro do sistema linear, retornara o vetor de residuos se o problema convergir
 * tol <- tolerancia para a norma do residuo
 * numiter <- retorna o numero de iteracoes ate a solucao ser atingida
 * Obs.: algoritmo pag 200 do livro do Yousef Saad
 * Link: http://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf
 */
bool solve_gradiente_conjugado(int N, double* A, double* x, double* b, double tol, int* numiter) {
	double* p = (double*) malloc(N * sizeof(double));
	double* q = (double*) malloc(N * sizeof(double));
	//r <- b - A*x
	cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, -1.0, A, N, x, 1, 1.0, b, 1);

	//p <- r
	cblas_dcopy(N, b, 1, p, 1);

	for (int j = 0; j < N; j++)
		q[j] = 0;

	double alpha1 = cblas_ddot(N, b, 1, b, 1);
	int convergiu = FALSE;
	int i;
	for (i = 0; (!convergiu) && (i < MAXITER); i++) {
		//q_j <- A*p_j
		cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, A, N, p, 1, 0.0, q, 1);

		//alpha_j = <r_j, r_j>/<A*p_j, p_j>
		double alpha2 = cblas_ddot(N, q, 1, p, 1);
		double alpha = alpha1 / alpha2;

		//x_{j+1} <- x_j + alpha_j*p_j
		cblas_daxpy(N, alpha, p, 1, x, 1);

		//r_{j+1} <- r_j - alpha_j*A*p_j
		cblas_daxpy(N, -alpha, q, 1, b, 1);

		double res = cblas_dnrm2(N, b, 1);
		if (res < tol)
			convergiu = TRUE;
		else {
			//beta_j = <r_{j+1}, r_{j+1}>/<r_j}, r_j>
			double beta1 = cblas_ddot(N, b, 1, b, 1);
			double beta2 = alpha1;
			double beta = beta1 / beta2;

			//p_{j+1} <- r_{j+1} + beta_j*p_j
			cblas_daxpy(N, 1.0 / beta, b, 1, p, 1);
			cblas_dscal(N, beta, p, 1);

			alpha1 = beta1;
		}
	}

	free(q);
	free(p);

	*numiter = i;

	return convergiu;
}

void reinicializa_vetores(int N, double* b, double* b0, double* x, double* x0) {
	memcpy(b, b0, sizeof(double) * N);
	memcpy(x, x0, sizeof(double) * N);
}

int main() {
	double tol = 1.0e-6;
	int N = 3;
	double A[9] = { 2, -1, 0, -1, 2, -1, 0, -1, 2 };
	double b0[3] = { 5, -7, 7 };
	double x0[3] = { -1, -1, -1 };

	double b[3];
	double x[3];
	int numiter;
	bool convergiu;

	printf("Steepest descent\n");
	reinicializa_vetores(N, b, b0, x, x0);
	convergiu = solve_steepest_descent(N, A, x, b, tol, &numiter);
	print_solucao(convergiu, numiter, N, x);

	printf("Gradiente conjugado\n");
	reinicializa_vetores(N, b, b0, x, x0);
	convergiu = solve_gradiente_conjugado(N, A, x, b, tol, &numiter);
	print_solucao(convergiu, numiter, N, x);

	return 0;
}
