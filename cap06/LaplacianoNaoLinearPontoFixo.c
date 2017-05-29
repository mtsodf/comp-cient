#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <math.h>
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
		printf("%10.4lf", v[i]);
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
double u(double x, double y) {
	//return exp(x) + exp(y);
	//return sin(M_PI * x) * sin(M_PI * y);
	return x + y;
}

double ld(double x, double y) {
	//return exp(x) + exp(y);
	//return -2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
	return 4*(x+y);
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

void printMatriz(matriz A) {
	for (int i = 0; i < A.nd; ++i) {
		printf("Diagonal %d = ", A.desloc[i]);
		printVec(A.n - A.desloc[i], A.val[i]);
	}
}

double trans(double T){
	return 1+T*T;
	//return 1.0;
}

int m2g(int i, int j, int n){
	return j*(n-1) + i;
}


void free_matriz(matriz *A){
	for (int i = 0; i < A->nd; ++i)
	{
		free(A->val[i]);
	}

	free(A->desloc);
}


void laplace(int n, double * x, double* xsol) {

	double h = 1.0 / n;

	double rnorm;

	int unknows = (n - 1) * (n - 1);

	FILE* saida;

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

	cblas_dcopy (unknows, x, 1, xsol, 1);


	b = (double*) malloc(sizeof(double) * (n - 1) * (n - 1));

	for (int i = 0; i < (unknows); ++i){
		A.val[0][i] = -4 / (h * h);
	}

	for (int i = 0; i < (unknows - 1); ++i) {
		if (i % (n - 1) == n - 2) {
			A.val[1][i] = 0.0;
		} else {
			A.val[1][i] = 1 / (h * h);
		}

	}
	for (int i = 0; i < (unknows - n + 1); ++i){
		A.val[2][i] = 1 / (h * h);
	}

	// Criando lado direito
	for (int j = 0; j < n - 1; ++j) {
		for (int i = 0; i < n - 1; ++i) {
			double x = (i + 1) * h;
			double y = (j + 1) * h;

			int ind = m2g(i,j,n);

			b[ind] = ld(x, y)/ trans(xsol[ind]);

			//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		}
	}

	// Bordas x=0 e x=1
	for (int j = 0; j < n - 1; ++j) {
		int ind = m2g(0,j,n);
		double x = 0.0;
		double y = (j + 1) * h;

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= u(x, y) / (h * h);

		ind = m2g(n-2, j, n);
		x = 1.0;

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= u(x, y) / (h * h);

	}

	// Bordas y=0 e y=1
	for (int i = 0; i < n - 1; ++i) {
		int ind = i;
		double x = (i + 1) * h;
		double y = 0.0;

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= u(x, y) / (h * h);

		ind = m2g(i,n-2,n);
		y = 1.0;

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= u(x, y) / (h * h);

	}

	//Parcela Relativa a d(1+T*T) * dT
	for (int j = 0; j < n - 1; ++j) {
		for (int i = 0; i < n - 1; ++i) {
			double x = (i + 1) * h;
			double y = (j + 1) * h;

			double T1, T2;

			int ind = m2g(i,j,n);

			if(i == 0){
				T1 = u(0.0, y);
			} else{
				T1 = xsol[m2g(i-1,j,n)];
			}

			if(i==n-2){
				T2 = u(1.0,y);
			}else{
				T2 = xsol[m2g(i+1,j,n)];
			}

			double aux = (trans(T2)-trans(T1))/(2*h);
			aux *= (T2 - T1)/(2*h);

			//aux = 2*xsol[ind] ;
			//aux *= (T2 - T1)/(2*h);
			//aux *= (T2 - T1)/(2*h);
			b[ind] -= aux/trans(xsol[ind]);

			if(j == 0){
				T1 = u(x, 0.0);
			} else{
				T1 = xsol[m2g(i,j-1,n)];
			}

			if(j==n-2){
				T2 = u(x,1.0);
			}else{
				T2 = xsol[m2g(i,j+1,n)];
			}
			
			aux = (trans(T2)-trans(T1))/(2*h);
			aux *= (T2 - T1)/(2*h);
			//aux = 2*xsol[ind] ;
			//aux *= (T2 - T1)/(2*h);
			//aux *= (T2 - T1)/(2*h);
			b[ind] -= aux/trans(xsol[ind]);
		}
	}	


	cg(A, b, xsol, &rnorm);
	
	saida = fopen("saida.txt", "w");

	fprintf(saida, "n = %d\n", n);

	fprintf(saida, "\nSolucao Encontrada\n");
	for (int i = 0; i < n - 1; ++i) {
		printVecFile(n - 1, xsol + i * (n - 1), saida);
	}


	double * sol = (double *) malloc(sizeof(double) * unknows);
	calc_sol(n, sol);

	fprintf(saida, "\nSolucao Real\n");
	for (int i = 0; i < n - 1; ++i) {
		printVecFile(n - 1, sol + i * (n - 1), saida);
	}

	cblas_daxpy(unknows, -1.0, xsol, 1, sol, 1);

	double difsol = cblas_dnrm2(unknows, sol, 1);

	fprintf(saida, "\n\tNorma Residuo = %e\n", rnorm);
	printf("\n\tNorma Residuo = %e\n", rnorm);
	printf("\n\tNorma Solucao = %e\n", difsol);
	fprintf(saida, "\n\tNorma Solucao = %e\n", difsol);

	fclose(saida);

	free_matriz(&A);
	free(sol);

}



int main(int argc, char **argv) {

	int n = atoi(argv[1]);
	printf("%s\n", argv[1]);
	printf("n = %d \n", n);
	

	int unknows = (n-1)*(n-1);

	double *x1 = (double*)malloc(sizeof(double)*unknows);
	double *x2 = (double*)malloc(sizeof(double)*unknows);
	double *aux;


	for (int i = 0; i < unknows; ++i)
	{
		x1[i] = 1.0;	
	}

	double rnorm = 1.0;
	for (int i = 0; i < MAXITER && rnorm > 1e-10; ++i)
	{
		laplace(n, x1, x2);

		cblas_daxpy (unknows, -1.0, x2, 1, x1, 1);

		rnorm = cblas_dnrm2 (unknows, x1, 1);

		printf("\t%d rnorm  = %lf\n", i, rnorm);

		aux = x2;
		x2 = x1;
		x1 = aux;
	}
	
	free(x1);
	free(x2);

}
