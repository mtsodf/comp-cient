/*
 * Laplaciano.c
 *
 *  Created on: 26 de mai de 2017
 *      Author: cx3d
 */

#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <math.h>
#include <string.h>
#include "solvers.h"

/**
 * Comentários:
 *
 * - steppest descent tem convergencia mais lenta (maior numero de iteracoes)
 * - se inserir as condicoes de contorno na matriz, ela deixa de ser simetrica
 * - precondicionamento nao fez diferenca (no problema linear)
 */

void printMatriz(matriz);

void printVec(int n, double *v) {
	for (int i = 0; i < n; ++i) {
		printf("%10.4e ", v[i]);
	}

	printf("\n");
}

void printVecFile(int n, double *v, FILE* f) {
	for (int i = 0; i < n; ++i) {
		fprintf(f, "%10.4lf", v[i]);
	}

	fprintf(f, "\n");
}

/**
 * Solucao analitica
 */
double u_sol(double x, double y, double t) {
	return 0.0;
	return 20*x * (1 - x) * y * (1 - y) * exp(-t);

}

double contorno(double x, double y, double t) {
	if(y>=1){
		return 1.0;
	}
	return 0.0;
}

double ld(double x, double y, double t) {
	return 0.0;
	double aux = -2 * (1 - x) * x * exp(-t) - 2 * (1 - y) * y * exp(-t)
			 + (1 - x) * x * (1 - y) * y * exp(-t);
	return 20*aux;
}

void calc_sol(int n, double *sol, double t) {
	double h = 1.0 / n;
	int ind = 0;
	double x, y;

	for (int j = 0; j < n - 1; ++j) {
		y = (j + 1) * h;
		for (int i = 0; i < n - 1; ++i) {
			x = (i + 1) * h;
			sol[ind] = u_sol(x, y, t);

			ind++;
		}
	}
}

void metodo_implicito(double t, double dt, double* u0, int n) {

	double h = 1.0 / n;

	double rnorm;

	int unknows = (n - 1) * (n - 1);

	matriz A;

	//Alocando Matriz do Laplaciano
	A.val = (double**) malloc(sizeof(double*) * 3);
	A.val[0] = (double *) malloc(sizeof(double) * unknows);
	A.val[1] = (double *) malloc(sizeof(double) * (unknows - 1));
	A.val[2] = (double *) malloc(sizeof(double) * (unknows - n + 1));
	A.desloc = (int *) malloc(sizeof(int) * (3));
	A.desloc[0] = 0;
	A.desloc[1] = 1;
	A.desloc[2] = n - 1;
	A.nd = 3;
	A.n = unknows;

	double * b;

	t += dt;

	b = (double*) malloc(sizeof(double) * unknows);

	for (int i = 0; i < (unknows); ++i)
		A.val[0][i] = -4 / (h * h) - 1/dt;
	for (int i = 0; i < (unknows - 1); ++i) {
		if (i % (n - 1) == n - 2) {
			A.val[1][i] = 0.0;
		} else {
			A.val[1][i] = 1 / (h * h);
		}

	}

	for (int i = 0; i < (unknows - n + 1); ++i)
		A.val[2][i] = 1 / (h * h);

	// Criando lado direito
	for (int j = 0; j < n - 1; ++j) {
		for (int i = 0; i < n - 1; ++i) {
			double x = (i + 1) * h;
			double y = (j + 1) * h;

			int ind = j * (n - 1) + i;

			b[ind] = ld(x, y, t) - u0[ind]/dt;
		}
	}

	// Bordas x=0 e x=1
	for (int j = 0; j < n - 1; ++j) {
		int ind = j * (n - 1) + 0;
		double x = 0.0;
		double y = (j + 1) * h;

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= contorno(x, y, t) / (h * h);

		ind = j * (n - 1) + n - 2;
		x = 1.0;

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= contorno(x, y, t) / (h * h);

	}

	// Bordas y=0 e y=1
	for (int i = 0; i < n - 1; ++i) {
		int ind = i;
		double x = (i + 1) * h;
		double y = 0.0;

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= contorno(x, y, t) / (h * h);

		ind = (n - 2) * (n - 1) + i;
		y = 1.0;

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= contorno(x, y, t) / (h * h);

	}


	cg(A, b, u0, &rnorm);

	printf("\tNorma CG = %e\n", rnorm);

}

void define_vizinhos(double t, int idx, double* T, int n, double* T_O,
		double* T_N, double* T_S, double* T_E, double* T_W) {
	double h = 1.0 / n;
	int i = idx % (n - 1);
	int j = idx / (n - 1);
	double x = (i + 1) * h;
	double y = (j + 1) * h;

	int idx_N = i + (j + 1) * (n - 1);
	int idx_S = i + (j - 1) * (n - 1);
	int idx_E = i + 1 + j * (n - 1);
	int idx_W = i - 1 + j * (n - 1);

	*T_O = T[idx];

	if (j == n - 2)
		*T_N = contorno(x, 1, t);
	else
		*T_N = T[idx_N];

	if (j == 0)
		*T_S = contorno(x, 0, t);
	else
		*T_S = T[idx_S];

	if (i == n-2)
		*T_E = contorno(1, y, t);
	else
		*T_E = T[idx_E];

	if (i == 0)
		*T_W = contorno(0, y, t);
	else
		*T_W = T[idx_W];
}

void metodo_explicito(double t, double dt, double* u0, double* u, int n) {
	int unknowns = (n - 1) * (n - 1);
	double h = 1.0 / n;

	for (int ind = 0; ind < unknowns; ind++) {
		int i = ind % (n - 1);
		int j = ind / (n - 1);

		double x = (i + 1) * h;
		double y = (j + 1) * h;

		//laplaciano
		double u_C, u_N, u_S, u_E, u_W;
		define_vizinhos(t, ind, u0, n, &u_C, &u_N, &u_S, &u_E, &u_W);
		double lapl = (u_W - 2 * u_C + u_E) / (h * h)
				+ (u_N - 2 * u_C + u_S) / (h * h);
		u[ind] = (lapl - ld(x, y, t)) * dt + u_C;
	}
}

void comparar_sol_analitica(int unknowns, double* u, double* sol) {
	cblas_daxpy(unknowns, -1.0, u, 1, sol, 1);
	double norm = cblas_dnrm2(unknowns, sol, 1);
	printf("\tNorma de |u - sol| = %e \n", norm);
}

void printMatriz(matriz A) {
	for (int i = 0; i < A.nd; ++i) {
		printf("Diagonal %d = ", A.desloc[i]);
		printVec(A.n - A.desloc[i], A.val[i]);
	}
}

void printSolucao(FILE* saida, int n, double* sol){
	int N = n - 1;
	int ind = 0;
	for (int j = 0; j < n - 1; ++j) {
		for (int i = 0; i < n - 1; ++i) {
			fprintf(saida, "%f ", sol[ind++]);
		}
		fprintf(saida, "\n");
	}

}

int main(int argc, char **argv) {

	int n = atoi(argv[1]);
	printf("%s\n", argv[1]);
	printf("n = %d \n", n);

	int metodo = 0; double *aux;


	int unknowns = (n - 1) * (n - 1);

	double* u0 = (double*) malloc(sizeof(double) * unknowns);
	double* u = (double*) malloc(sizeof(double) * unknowns);
	double* sol = (double*) malloc(sizeof(double) * unknowns);

	//loop do metodo explicito
	double dt = 0.01;
	double t = 0;

	if(argc > 2){
		dt = atof(argv[2]);
	}

	int maxiters = 10;
	if(argc > 3){
		maxiters = atoi(argv[3]);
	}	

	//condicao inicial
	calc_sol(n, u0, 0);

	FILE* saida;

	saida = fopen("saida.txt", "w");

	fprintf(saida, "%d\n", n);


	for (int iter = 0; iter < maxiters; iter++) {

		printf("t = %lf\n", t+dt);
		fprintf(saida, "%lf\n", t+dt);
		//metodo_implicito(t, dt, u0, n);

		metodo_explicito(t, dt, u0, u, n);
		aux = u0;
		u0 = u;
		u = aux;

		t += dt;
		calc_sol(n, u, t);
		comparar_sol_analitica(unknowns, u0, u);
		printSolucao(saida, n, u0);

	}

	fclose(saida);

	//laplace(n);
}

