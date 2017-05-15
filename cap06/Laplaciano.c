#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <math.h>
#include "solvers.h"




void printMatriz(matriz);

void printVec(int n, double *v){
	for (int i = 0; i < n; ++i)
	{
		printf("%10.4lf", v[i]);
	}

	printf("\n");
}

void printVecFile(int n, double *v, FILE* f){
	for (int i = 0; i < n; ++i)
	{
		fprintf(f, "%10.4lf", v[i]);
	}

	fprintf(f, "\n");
}



double ld(double x, double y){
	return exp(x) + exp(y);
	//return 0.0;
	//return exp(x);
}



void set_element(matriz A, int i, int j, double v){
	int aux;

	if(j < i){
		aux = i;
		i = j;
		j = aux;
	}

	int desloc = j - i;

	for (int k = 0; k < A.nd; ++k)
	{
		if(A.desloc[k] == desloc){
			A.val[k][i] = v; 
			return ;
		}
	}

	printf("Setando valor fora das diagonais de A\n");
	v = 1.0 / 0.0;
}





void calc_sol(int n, double *sol){
	double h = 1.0 / n;
	int ind = 0;
	double x, y;

	for (int j = 0; j < n - 1; ++j)
	{	
		y = (j+1)*h;
		for (int i = 0; i < n - 1 ; ++i)
		{
			x = (i+1)*h;
			sol[ind] = exp(x) + exp(y);

			ind++;
		}
	}
}

void printMatriz(matriz A){
	for (int i = 0; i < A.nd; ++i)
	{	
		printf("Diagonal %d = ", A.desloc[i]);
		printVec(A.n-A.desloc[i], A.val[i]);
	}
}

void laplace(int n){

	double * x = (double *) malloc(sizeof(double) * (n-1)*(n-1));

	double h = 1.0 / n;

	double rnorm;

	int unknows = (n-1)*(n-1);

	FILE* saida;

	matriz A;

	//Alocando Matriz do Laplaciano
	A.val = (double**) malloc(sizeof(double*)*3);
	A.val[0] = (double *) malloc(sizeof(double)*unknows);
	A.val[1] = (double *) malloc(sizeof(double)*(unknows-1));
	A.val[2] = (double *) malloc(sizeof(double)*(unknows-n+1));
	A.desloc = (int *) malloc(sizeof(int)*(3));
	A.desloc[0] = 0;
	A.desloc[1] = 1;
	A.desloc[2] = n-1;
	A.nd = 3;
	A.n = unknows;

	double * b;

	b = (double*) malloc(sizeof(double) * (n-1)*(n-1));

	for (int i = 0; i < (unknows); ++i) A.val[0][i] = -4/(h*h) ;


	for (int i = 0; i < (unknows-1); ++i) {
		if(i%(n-1) == n-2){
			A.val[1][i] = 0.0 ;
		} else{
			A.val[1][i] = 1/(h*h) ;
		}
		
	}
	for (int i = 0; i < (unknows-n+1); ++i) A.val[2][i] = 1/(h*h) ;


	// Criando lado direito
	for (int j = 0; j < n - 1; ++j)
	{
		for (int i = 0; i < n - 1; ++i)		
		{	
			double x = (i+1) * h;
			double y = (j+1) * h;

			int ind = j*(n-1) + i;

			b[ind] = ld(x,y);

			//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		}
	}

	// Bordas x=0 e x=1
	for (int j = 0; j < n - 1; ++j)
	{
		int ind = j*(n-1) + 0;	
		double x = 0.0;
		double y = (j+1) * h;	

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= ld(x, y)/(h*h);

		ind = j*(n-1) + n - 2;	
		x = 1.0;

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= ld(x,y)/(h*h);

	}	

	// Bordas y=0 e y=1
	for (int i = 0; i < n - 1; ++i)
	{
		int ind = i ;	
		double x = (i+1) * h;
		double y = 0.0;	

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= ld(x, y)/(h*h);

		ind = (n-2)*(n-1) + i;
		y = 1.0;

		//printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= ld(x,y)/(h*h);

	}		

	double *x0 = (double*) malloc(sizeof(double)*A.n);
	for (int i = 0; i < A.n; ++i)
	{
		x0[i] = 1.0;
	}
	
	//cg_precon(A, b, x0, &rnorm);
	int niters;
	solve_steepest_descent(A.n, A, x0, b, 1e-10, &niters);
	
	saida = fopen("saida.txt", "w");

	fprintf(saida, "n = %d\n", n);

	fprintf(saida, "\nSolucao Encontrada\n");
	for (int i = 0; i < n - 1; ++i)
	{
		printVecFile(n-1, x0 + i*(n-1), saida);
	}
	
	double * sol = (double *) malloc(sizeof(double)*unknows);
	calc_sol(n, sol);


	fprintf(saida, "\nSolucao Real\n");
	for (int i = 0; i < n - 1; ++i)
	{
		printVecFile(n-1, sol + i*(n-1), saida);
	}

	cblas_daxpy (unknows, -1.0, x0, 1, sol, 1);

	fprintf(saida, "\n\tNorma Residuo = %lf\n", rnorm);
	printf("\n\tNorma Residuo = %lf\n", rnorm);

	fclose(saida);

}



int main(int argc, char **argv){
	
	int n = atoi(argv[1]);
	printf("%s\n", argv[1]);
	printf("n = %d \n", n);
	laplace(n);

}