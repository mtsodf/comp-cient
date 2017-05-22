#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

#define dT  0.0001

/**
 * Comentários:
 *
 * - steppest descent tem convergencia mais lenta (maior numero de iteracoes)
 * - se inserir as condicoes de contorno na matriz, ela deixa de ser simetrica
 * - precondicionamento nao fez diferenca (no problema linear)
 */

/**
 * Solucao analitica
 */
double u(double x, double y) {
	//return exp(x) + exp(y);
	return sin(M_PI * x) * sin(M_PI * y);
}

double ld(double x, double y) {
	//return exp(x) + exp(y);
	return -2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
}

void calc_sol(int n, double *sol) {
	double h = 1.0 / n;
	int ind = 0;
	double x, y;

	for (int j = 0; j < n - 1; ++j) {
		y = (j + 1) * h;
		for (int i = 0; i < n - 1; ++i) {
			x = (i + 1) * h;
			sol[ind] = u(x, y);

			ind++;
		}
	}
}

/*
 * Calcula o valor de k na interface
 */
double k_interface(double T1, double T2) {
	//return ((1 + T1 * T1) + (1 + T2 * T2)) / 2;
	return 1.0;
}

double contorno(double x, double y) {
	return 0.0;
}

double residuo(int n, gsl_vector *T, int i, int j) {

	double r;

	double Tij, Tv;

	int ind, indv; //indice do vizinho

	double x, y;

	double h = 1.0 / n; //distancia dos pontos

	ind = i + (n - 1) * j;
	Tij = gsl_vector_get(T, ind);

	r = 0;

	if (i > 0) {
		indv = i - 1 + (n - 1) * j;
		Tv = gsl_vector_get(T, indv);
	} else {
		x = 0.0;
		y = (j + 1) * h;
		Tv = contorno(x, y);
	}

	r += k_interface(Tv, Tij) * (Tv - Tij);

	if (i < n - 2) {
		indv = i + 1 + (n - 1) * j;
		Tv = gsl_vector_get(T, indv);
	} else {
		x = 1.0;
		y = (j + 1) * h;
		Tv = contorno(x, y);
	}

	r += k_interface(Tv, Tij) * (Tv - Tij);

	if (j > 0) {
		indv = i + (n - 1) * (j - 1);
		Tv = gsl_vector_get(T, indv);
	} else {
		x = (i + 1) * h;
		y = 0.0;
		Tv = contorno(x, y);
	}

	r += k_interface(Tv, Tij) * (Tv - Tij);

	if (j < n - 2) {
		indv = i + (n - 1) * (j + 1);
		Tv = gsl_vector_get(T, indv);
	} else {
		x = (i + 1) * h;
		y = 1.0;
		Tv = contorno(x, y);
	}

	r += k_interface(Tv, Tij) * (Tv - Tij);

	r /= h * h;

	return r;
}

double d_residuo(int n, gsl_vector *T, int i, int j, int ind) {
	double r1, r2;

	r1 = residuo(n, T, i, j);

	//T[ind] += dT;
	gsl_vector_set(T, ind, gsl_vector_get(T, ind) + dT);
	r2 = residuo(n, T, i, j);

	//T[ind] -= dT;
	gsl_vector_set(T, ind, gsl_vector_get(T, ind) - dT);

	return (r2 - r1) / dT;

}

void calc_residuo(int n, gsl_vector *T, double *r) {
	int unknows = (n - 1) * (n - 1); //numero de variaveis

	double x, y;

	int ind = 0;

	for (int j = 0; j < n - 1; ++j)
		for (int i = 0; i < n - 1; i++) {
			{
				ind = i + (n - 1) * j;
				r[ind] = residuo(n, T, i, j);
			}
		}

}

void calc_jacobiano(gsl_spmatrix* A, int n, gsl_vector *T) {

	int unknows = (n - 1) * (n - 1);

	int i, j;

	for (int linha = 0; linha < unknows; linha++) {
		i = linha % (n - 1);
		j = linha / (n - 1);
		int linha_N = i + (j - 1) * (n - 1);
		int linha_S = i + (j + 1) * (n - 1);
		int linha_E = i + 1 + j * (n - 1);
		int linha_W = i - 1 + j * (n - 1);
		gsl_spmatrix_set(A, linha, linha, d_residuo(n, T, i, j, linha));
		if (i > 0)
			gsl_spmatrix_set(A, linha, linha_W, d_residuo(n, T, i, j, linha_W));
		if (i < n - 2)
			gsl_spmatrix_set(A, linha, linha_E, d_residuo(n, T, i, j, linha_E));
		if (j > 0)
			gsl_spmatrix_set(A, linha, linha_N, d_residuo(n, T, i, j, linha_N));
		if (j < n - 2)
			gsl_spmatrix_set(A, linha, linha_S, d_residuo(n, T, i, j, linha_S));
	}
}

int main(int argc, char **argv) {

	int n = atoi(argv[1]);
	printf("%s\n", argv[1]);
	printf("n = %d \n", n);

	int unknows = (n - 1) * (n - 1);
	gsl_vector* T = gsl_vector_alloc(unknows);
	for (int i = 0; i < unknows; ++i) {
		gsl_vector_set(T, i, i);
	}

	gsl_spmatrix* A = gsl_spmatrix_alloc_nzmax(unknows, unknows, 5 * unknows,
	GSL_SPMATRIX_TRIPLET);
	calc_jacobiano(A, n, T);

	gsl_spmatrix_fprintf(stdout, A, "%5.3f");

}
