#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <math.h>

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
		printf("%6.2lf", v[i]);
	}

	printf("\n");
}

double ld(double x, double y){
	return exp(-x) + exp(-y);
	//return 0.0;
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

	for (int i = 0; i < (unknows); ++i) A.val[0][i] = -4 /(h*h) ;
	for (int i = 0; i < (unknows-1); ++i) A.val[1][i] = -1 /(h*h) ;
	for (int i = 0; i < (unknows-n+1); ++i) A.val[2][i] = -1 /(h*h) ;





	for (int j = 0; j < n - 1; ++j)
	{
		for (int i = 0; i < n - 1; ++i)		
		{	
			double x = (i+1) * h;
			double y = (j+1) * h;

			int ind = j*(n-1) + i;

			b[ind] = ld(x,y);
		}
	}

	for (int j = 0; j < n - 1; ++j)
	{
		int ind = j*(n-1);	
		double x = 0.0;
		double y = (j+1) * h;	

		b[ind] -= ld(x, y);

		ind = j*(n-1) + n - 2;	
		x = 1.0;
		b[ind] -= ld(x,y);

	}	

	for (int i = 0; i < n - 1; ++i)
	{
		int ind = i ;	
		double x = (i+1) * h;
		double y = 0.0;	

		b[ind] -= ld(x, y);

		ind = (n-2)*(n-1) + i;	
		y = 1.0;
		b[ind] -= ld(x,y);

	}		

	printf("b = ");
	printVec(unknows, b);
}


void matmult(int n, matriz A, double *b, double *r){

	double *aux = (double *) malloc(sizeof(double)*n);




	for (int i = 0; i < A.nd; ++i)
	{	
		int dsl = A.desloc[i];
		vdMul(n-dsl, A.val[i], b+dsl, aux);
		printf("aux%d = ", i);
		printVec(n-dsl, aux);
		if(i > 0){
			//Diagonal Superior
			cblas_daxpy (n-dsl, 1.0, aux, 1, r, 1);

			//Diagonal Inferior
			vdMul(n-dsl, A.val[i], b, aux);
			cblas_daxpy (n-dsl, 1.0, aux, 1, r+dsl, 1);

		} else{
			cblas_dcopy (n, aux, 1, r, 1);
		}



		printf("r%d = ", i);
		printVec(n, r);
		
	}

	printf("r = ");
	printVec(n, r);

}


int main(){
	

	matriz A;
	double b[4], r[4];
	int n, nd;
	int aux;

	


	FILE *fp;

	fp = fopen("../entrada.txt", "r");
	fscanf(fp, "%d %d", &n, &nd);
	printf("n = %d\tnd = %d\n", n, nd);

	A.val = (double **) malloc(sizeof(double*) * nd);
	A.desloc = (int*) malloc(sizeof(int) * nd);
	A.n = n;
	A.nd = nd;

	for (int i = 0; i < nd; ++i)
	{
		fscanf(fp, "%d", &A.desloc[i]);
		A.val[i] = (double *) malloc(sizeof(double) * (n - A.desloc[i]));
		for (int j = 0; j < n - A.desloc[i]; ++j)
		{
			fscanf(fp, "%lf", &A.val[i][j]);
		}
	}

	for (int i = 0; i < A.nd; ++i)
	{	
		printf("Diagonal %d = ", A.desloc[i]);
		printVec(n-A.desloc[i], A.val[i]);
	}
	

	for (int i = 0; i < n; ++i)   fscanf(fp, "%lf", &b[i]);

	matmult(n, A, b, r);

	printVec(n, r);

	for (int i = 0; i < A.nd; ++i) free(A.val[i]);
	

	free(A.val);
	free(A.desloc);


	fclose(fp);

	laplace(3);

	//vdMul( 2, a+1, b, y );

}