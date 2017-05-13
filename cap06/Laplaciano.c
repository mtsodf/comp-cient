#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <math.h>


/* Estrutura para guardar as matrizes 
 *
 *
*/
struct struct_matriz 
{
	double** val;
	int* desloc;
	int n;
	int nd;
};


typedef struct struct_matriz matriz;


void printVec(int n, double *v){
	for (int i = 0; i < n; ++i)
	{
		printf("%10.4lf", v[i]);
	}

	printf("\n");
}


double ld(double x, double y){
	return exp(x) + exp(y);
	//return 0.0;
	//return exp(x);
}

/* Calcula o produto A*b e retorna em r.
*  O vetor r deve estar previamente alocado.
*/
void matmult(int n, matriz A, double *b, double *r){

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


void cg(matriz A, double* b, double * x0){
	double *r, *p;
	double *Ap;
	double rnorm, alpha, rdotr, rdotr2, beta;

	Ap = (double*) malloc(sizeof(double)*A.n);
	 p = (double*) malloc(sizeof(double)*A.n);
	r  = b;


	matmult(A.n, A, x0, Ap);

	cblas_daxpy (A.n, -1.0, Ap, 1, r, 1);

	rnorm = cblas_dnrm2 (A.n, r, 1);
	cblas_dcopy (A.n, r, 1, p, 1);

	while(rnorm > 1e-10){
		matmult(A.n, A, p, Ap);

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

		rnorm = sqrt(rdotr2);
	}

	free(Ap);
	
}

void laplace(int n){

	double * x = (double *) malloc(sizeof(double) * (n-1)*(n-1));

	double h = 1.0 / n;

	int unknows = (n-1)*(n-1);

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

	for (int i = 0; i < (unknows); ++i) A.val[0][i] = -4 /(2*h) ;
	for (int i = 0; i < (unknows-1); ++i) {
		if(i%(n-1) == n- 2){
			A.val[1][i] = 0.0 ;
		} else{
			A.val[1][i] = 1/(2*h) ;
		}
		
	}
	for (int i = 0; i < (unknows-n+1); ++i) A.val[2][i] = 1/(2*h) ;


	// Criando lado direito
	for (int j = 0; j < n - 1; ++j)
	{
		for (int i = 0; i < n - 1; ++i)		
		{	
			double x = (i+1) * h;
			double y = (j+1) * h;

			int ind = j*(n-1) + i;

			b[ind] = ld(x,y);

			printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		}
	}

	// Bordas x=0 e x=1
	for (int j = 0; j < n - 1; ++j)
	{
		int ind = j*(n-1) + 0;	
		double x = 0.0;
		double y = (j+1) * h;	

		printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= ld(x, y)/(2*h);

		ind = j*(n-1) + n - 2;	
		x = 1.0;

		printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= ld(x,y)/(2*h);

	}	

	// Bordas y=0 e y=1
	for (int i = 0; i < n - 1; ++i)
	{
		int ind = i ;	
		double x = (i+1) * h;
		double y = 0.0;	

		printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= ld(x, y)/(2*h);

		ind = (n-2)*(n-1) + i;	
		y = 1.0;

		printf("x=%lf\ty=%lf\tf=%lf\n", x, y, ld(x,y));
		b[ind] -= ld(x,y)/(2*h);

	}		

	printf("b = ");
	printVec(unknows, b);

	for (int i = 0; i < A.nd; ++i)
	{	
		printf("Diagonal %d = ", A.desloc[i]);
		printVec(A.n-A.desloc[i], A.val[i]);
	}

	double *x0 = (double*) malloc(sizeof(double)*A.n);
	for (int i = 0; i < A.n; ++i)
	{
		x0[i] = 1.0;
	}
	
	cg(A, b, x0);
	
	for (int i = 0; i < n - 1; ++i)
	{
		printVec(n-1, x0 + i*(n-1));
	}
	
}



int main(int argc, char **argv){
	
	int n = atoi(argv[1]);
	printf("%s\n", argv[1]);
	printf("n = %d \n", n);
	laplace(n);

}