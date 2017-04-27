//============================================================================
// Name        : RayleightCoeffIter.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

using namespace std;

void alocar(double ** x, int n){
	*x = (double*) malloc(sizeof(double)*n);
}

void printVec(double *x, int n){
	for (int i = 0; i < n; ++i) {
		printf("%8.4lf", x[i]);
	}
	//printf("\n");
}

void printMat(double *x, int n){
	int inic = 0;
	for (int i = 0; i < n; ++i) {
		printf("\t\t");
		printVec(&x[inic], n);
		printf("\n");
		inic += n;
	}

}

void RaileightCoefIter(double *A, double *x, int n, double *autovalor, int log){
	int info;

	double *y, *xant;

	double *Ac;

	double sigma, norm, aux;

	int *ipiv, nrhs = 1, passo = 0;


	alocar(&xant, n);
	ipiv = (int*) malloc(sizeof(int)*n);
	Ac = (double *)mkl_malloc( n*n*sizeof( double ), 64);

	for (int i = 0; i < n*n; ++i) {
		Ac[i] = A[i];
	}

	alocar(&y, n);

	do{


		cblas_dcopy (n, x, 1, xant, 1);

		cblas_dgemv (CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, x, 1, 0.0, y, 1);

		sigma = cblas_ddot (n, x, 1, y, 1);
		sigma = sigma /cblas_ddot(n, x, 1, x, 1);

		if(log){
			printf("%8d |", passo);
			printVec(x, n);
			printf(" |%8.4f\n", sigma);
		}


		for(int i = 0; i < n*n; ++i){
			Ac[i] = A[i];
		}

		int d = 0;
		for (int i = 0; i < n; ++i) {
			Ac[d] = A[d] - sigma;
			d += n + 1;
		}

		dgesv( &n, &nrhs, Ac, &n, ipiv, x, &n, &info );

        if( info > 0 ) {
        		cblas_dcopy (n, xant, 1, x, 1);
                break;
        }

		norm = x[cblas_idamax (n, x, 1)];
		norm = norm > 0 ? norm: -norm;

		for (int i = 0; i < n; ++i) {
			x[i] /= norm;
		}

		passo++;

		cblas_daxpy (n, -1.0, x, 1, xant, 1);

		norm = xant[cblas_idamax (n, xant, 1)];
		norm = norm > 0 ? norm: -norm;
	}while(norm > 1e-6);

	*autovalor = sigma;

	free(y);
	free(xant);
	mkl_free(Ac);

}

double randomDouble(){
	double random;
	random = rand();

	random -= RAND_MAX/2;

	random /= RAND_MAX/2;

	return random;
}

void RodarRandomico(double *A, int n, int iteracoes){
	double *x, autovalor;

	alocar(&x, n);

	for (int j = 0; j < iteracoes; ++j) {

		for (int i = 0; i < n; ++i) {
			x[i] = randomDouble();
		}

		RaileightCoefIter(A, x, n, &autovalor, 0);

		printVec(x, n); printf("   |%14.6lf", autovalor);
		printf("\n");

	}

	free(x);

}

int main() {
	int n;
	double *x, *A;
	double autovalor;
	FILE *fp;

	fp = fopen("entrada.txt", "r");
	fscanf(fp, "%d", &n);
	printf("n = %d\n", n);

	alocar(&x, n);

	A = (double *)mkl_malloc( n*n*sizeof( double ), 64);


	for (int i = 0; i < n*n; ++i) {
		fscanf(fp, "%lf", &A[i]);
	}

	for (int i = 0; i < n; ++i) {
		fscanf(fp, "%lf",&x[i]);
	}


	RodarRandomico(A, n, 20);

	return 0;

	printf("Inicializando Algoritmo com:\n");
	printf("\tMatriz A\n");
	printMat(A, n);
	printf("\tChute Inicial\n\t\t");
	printVec(x, n);printf("\n");


	printf("\n**********************************************\n");
	printf("              Inicio do Algoritmo\n");
	printf("**********************************************\n\n");
	RaileightCoefIter(A, x, n, &autovalor, 1);





	free(x);
	mkl_free(A);

	return 0;
}
