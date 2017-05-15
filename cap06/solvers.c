#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <math.h>
#include "solvers.h"



typedef struct struct_matriz matriz;


extern void printVec(int, double );

/* Calcula o produto A*b e retorna em r.
*  O vetor r deve estar previamente alocado.
*/
void mat_mult(int n, matriz A, double *b, double *r){

	double *aux = (double *) malloc(sizeof(double)*n);


	for (int i = 0; i < A.nd; ++i)
	{	
		int dsl = A.desloc[i];
		vdMul(n-dsl, A.val[i], b+dsl, aux);

		if(i > 0){
			//Diagonal Superior
			cblas_daxpy (n-dsl, 1.0, aux, 1, r, 1);

			//Diagonal Inferior
			vdMul(n-dsl, A.val[i], b, aux);
			cblas_daxpy (n-dsl, 1.0, aux, 1, r+dsl, 1);

		} else{
			cblas_dcopy (n, aux, 1, r, 1);
		}

	}

	free(aux);
}

void cg_precon(matriz A, double* b, double * x0, double *rnorm){
	double *r, *p, *z, *invD;
	double *Ap;
	double  alpha, rdotr, rdotr2, beta;

	Ap = (double*) malloc(sizeof(double)*A.n);
	invD = (double*) malloc(sizeof(double)*A.n);
	 p = (double*) malloc(sizeof(double)*A.n);
	 z = (double*) malloc(sizeof(double)*A.n);

	r  = b;

	for (int i = 0; i < A.n; ++i)
	{
		invD[i] = 1/A.val[0][i];
	}


	mat_mult(A.n, A, x0, Ap);

	cblas_daxpy (A.n, -1.0, Ap, 1, r, 1);

	*rnorm = cblas_dnrm2 (A.n, r, 1);
	cblas_dcopy (A.n, r, 1, z, 1);

	vdMul(A.n, invD, z, z);

	cblas_dcopy (A.n, z, 1, p, 1);


	int iters = 0;

	while(*rnorm > 1e-10){
		mat_mult(A.n, A, p, Ap);

		rdotr = cblas_ddot(A.n, r, 1, z, 1);

		alpha =  rdotr/cblas_ddot(A.n, Ap, 1, p, 1);

		cblas_daxpy (A.n, alpha, p, 1, x0, 1);
		cblas_daxpy (A.n, -alpha, Ap, 1, r, 1);

		vdMul(A.n, invD, r, z);

		rdotr2 = cblas_ddot(A.n, r, 1, z, 1);

		beta = rdotr2 / rdotr;

		//pj = beta*pj
		cblas_dscal (A.n, beta, p, 1);
 
		//$p_(j+1) = r + beta*pj$
		cblas_daxpy (A.n, 1.0, z, 1, p, 1);

		*rnorm = cblas_dnrm2 (A.n, r, 1);

		iters++;
	}

	printf("Convergiu em %d iteracoes\n", iters);

	free(Ap);
	free(p);
	free(z);
	free(invD);
	
}

void cg(matriz A, double* b, double * x0, double *rnorm){
	double *r, *p;
	double *Ap;
	double  alpha, rdotr, rdotr2, beta;

	Ap = (double*) malloc(sizeof(double)*A.n);
	 p = (double*) malloc(sizeof(double)*A.n);
	r  = b;


	mat_mult(A.n, A, x0, Ap);

	cblas_daxpy (A.n, -1.0, Ap, 1, r, 1);

	*rnorm = cblas_dnrm2 (A.n, r, 1);
	cblas_dcopy (A.n, r, 1, p, 1);

	int iters = 0;

	while(*rnorm > 1e-10){
		mat_mult(A.n, A, p, Ap);

		rdotr = cblas_ddot(A.n, r, 1, r, 1);

		alpha =  rdotr/ cblas_ddot(A.n, Ap, 1, p, 1);

		cblas_daxpy (A.n, alpha, p, 1, x0, 1);
		cblas_daxpy (A.n, -alpha, Ap, 1, r, 1);

		rdotr2 = cblas_ddot(A.n, r, 1, r, 1);

		beta = cblas_ddot(A.n, r, 1, r, 1) / rdotr;

		//pj = beta*pj
		cblas_dscal (A.n, beta, p, 1);
 
		//$p_(j+1) = r + beta*pj$
		cblas_daxpy (A.n, 1.0, r, 1, p, 1);

		*rnorm = sqrt(rdotr2);

		iters++;
	}

	printf("Convergiu em %d iteracoes\n", iters);

	free(Ap);
	
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
int solve_steepest_descent(int N, matriz A, double* x, double* b, double tol, int* numiter, double* res) {
	double* p = (double*) malloc(N * sizeof(double));
	double* Ap = (double*) malloc(N * sizeof(double));
	//r <- b - A*x
	mat_mult(N, A, x, Ap);
	cblas_daxpy (N, -1.0, Ap, 1, b, 1);

	//p <- A*r
	mat_mult(N, A, b, p);

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

		*res = cblas_dnrm2(N, b, 1);

		printf("Residuo = %f\n", *res);

		if (*res < tol)
			convergiu = TRUE;
		else
			//p <- A*r
			mat_mult(N, A, b, p);

	}

	free(p);

	*numiter = i;

	return convergiu;
}
