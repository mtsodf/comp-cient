#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

#define dT  0.0001



/*
 * Calcula o valor de k na interface
 */
double k_interface(double T1, double T2) {
	return ((1 + T1 * T1) + (1 + T2 * T2)) / 2;
	//return 1.0;
}


/**
 * Solucao analitica
 */
double u(double x, double y) {
	//return exp(x) + exp(y);
	return x + y;
}

double ld(double x, double y) {
	//return exp(x) + exp(y);	
	return 4*(x+y);
}

double contorno(double x, double y) {
	return u(x,y);
}

void calc_sol(int n, gsl_vector *sol) {
	double h = 1.0 / n;
	int ind = 0;
	double x, y;

	for (int j = 0; j < n - 1; ++j) {
		y = (j + 1) * h;
		for (int i = 0; i < n - 1; ++i) {
			x = (i + 1) * h;
			gsl_vector_set(sol, ind, u(x, y));
			ind++;
		}
	}
}

double residuo(int n, gsl_vector *T, int i, int j) {

	double r;

	double Tij, Tv;

	int ind, indv; //indice do vizinho

	double x, y;

	double h = 1.0 / n; //distancia dos pontos

	x = (i+1)*h;
	y = (j+1)*h;

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

	r -= ld(x,y);

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

void calc_residuo(int n, gsl_vector *T, gsl_vector *r) {
	int unknows = (n - 1) * (n - 1); //numero de variaveis

	double x, y;

	int ind = 0;

	for (int j = 0; j < n - 1; ++j)
		for (int i = 0; i < n - 1; i++) {
			{
				ind = i + (n - 1) * j;
				gsl_vector_set(r, ind, residuo(n, T, i, j));
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

void printVecFile(int n, gsl_vector *v, FILE* f) {
	for (int i = 0; i < n; ++i) {
		fprintf(f, "%8.3lf", gsl_vector_get(v, i));
	}

	fprintf(f, "\n");
}



int main(int argc, char **argv) {

	int n = atoi(argv[1]);
	printf("%s\n", argv[1]);
	printf("n = %d \n", n);

	int unknows = (n - 1) * (n - 1);
	int status, iter, newton_steps = 0;
	const double tol = 1.0e-6; /* solution relative tolerance */
	const size_t max_iter = 100; /* maximum iterations */
	const gsl_splinalg_itersolve_type *Solver = gsl_splinalg_itersolve_gmres;
	double rnorm;
	gsl_splinalg_itersolve *gmres_solver = gsl_splinalg_itersolve_alloc(Solver, unknows, 0);



	gsl_vector* T = gsl_vector_alloc(unknows);
	for (int i = 0; i < unknows; ++i) {
		gsl_vector_set(T, i, 0);
	}

	gsl_vector* r = gsl_vector_alloc(unknows);

	gsl_spmatrix* A = gsl_spmatrix_alloc_nzmax(unknows, unknows, 5 * unknows, GSL_SPMATRIX_TRIPLET);

	gsl_vector* dx = gsl_vector_alloc(unknows);

	gsl_vector_set_zero (dx);

	calc_residuo(n, T, r);
	gsl_blas_dscal (-1.0, r);
	//gsl_vector_fprintf (stdout, r, "%.5lf");

	rnorm = gsl_blas_dnrm2(r);
	printf("Residuo Inicial = %lf\n\n", rnorm);
	
	do {
		printf("Passo de Newton %d\n", newton_steps++);
		iter = 0;
		do{
			calc_jacobiano(A, n, T);
			status = gsl_splinalg_itersolve_iterate(A, r, tol, dx, gmres_solver);
		}while(status == GSL_CONTINUE && ++iter < max_iter);

		rnorm = gsl_splinalg_itersolve_normr(gmres_solver);
		printf("\tIteracoes Lineares = %d\n", iter);

		gsl_blas_daxpy(1.0, dx, T);

		calc_residuo(n, T, r);
		gsl_blas_dscal (-1.0, r);
		rnorm = gsl_blas_dnrm2(r);

		printf("\trnorm = %lf\n\n", rnorm);
	} while(rnorm>1e-6);

	FILE* saida = fopen("saida.txt", "w");

	fprintf(saida, "n = %d\n", n);

	fprintf(saida, "\nSolucao Encontrada\n");
	for (int i = 0; i < n - 1; ++i) {
		for (int j = 0; j < n - 1; ++j){
			fprintf(saida, "%8.3f", gsl_vector_get (T, i+j*(n-1)));
		}
		fprintf(saida, "\n");
	}

	gsl_vector * sol = gsl_vector_alloc(unknows);
	calc_sol(n, sol);

	fprintf(saida, "\nSolucao Real\n");
	for (int i = 0; i < n - 1; ++i) {
		for (int j = 0; j < n - 1; ++j){
			fprintf(saida, "%8.3f", gsl_vector_get (sol, i+j*(n-1)));
		}
		fprintf(saida, "\n");
	}	

	gsl_blas_daxpy(-1.0, T, sol);

	double difsol = gsl_blas_dnrm2(sol);

	fprintf(saida, "\n\tNorma Residuo = %e\n", difsol);
	printf("\n\tNorma Residuo = %e\n", difsol);

	fclose(saida);

	gsl_vector_free(sol);

	gsl_vector_free(T);
	gsl_vector_free(r);	
	gsl_vector_free(dx);	
	gsl_splinalg_itersolve_free(gmres_solver);

}
