#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include <sys/time.h>

#define dT  0.0001


void contar_tempo(int comeco_fim, int num){
	static double inicio[100];
	static double acumulado[100];

	
}


/*
 * Calcula o valor de k na interface
 */
double k_interface(double T1, double T2) {
	return ((1 + T1 * T1) + (1 + T2 * T2)) / 2;
	//return 1.0;
}

/*
 * Calcula o valor da derivada de k na interface
 */
double dk_interface(double T1, double T2) {
	return T1  + T2;
	//return 0.0;
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

void get_jacobi(gsl_spmatrix* A , gsl_vector* jacobi){
	for(int i = 0; i <A->size1; i++){
		gsl_vector_set(jacobi, i, 1/gsl_spmatrix_get(A, i, i));
	}
}

void apply_jacobi_mat(gsl_spmatrix* A , gsl_vector* jacobi){
	int i; 
	for(int k = 0; k < A->nz; k++){
		i = A->i[k];
		A->data[k] *= gsl_vector_get(jacobi, i);
	}
}

void apply_jacobi_vec(gsl_vector* x , gsl_vector* jacobi){
	for(int i = 0; i < x->size; i++){
		x->data[i] *= gsl_vector_get(jacobi, i);
	}
}
/**
* Dada uma posicao no vetor solucao, recupera o valor desse ponto
* e os valores dos nos vizinhos desse ponto na malha. Nos pontos proximos
* a fronteira, recupera os valores correspondentes no contorno.
*/
void define_vizinhos(int idx, gsl_vector* T, int n, double* T_O, double* T_N,
		double* T_S, double* T_E, double* T_W) {
	double h = 1.0 / n;
	int i = idx % (n - 1);
	int j = idx / (n - 1);
	double x = (i + 1) * h;
	double y = (j + 1) * h;

	int idx_N = i + (j - 1) * (n - 1);
	int idx_S = i + (j + 1) * (n - 1);
	int idx_E = i + 1 + j * (n - 1);
	int idx_W = i - 1 + j * (n - 1);

	*T_O = gsl_vector_get(T, idx);

	if ((i == 0) & (j == 0)) {
		*T_N = contorno(x, 0);
		*T_S = gsl_vector_get(T, idx_S);
		*T_E = gsl_vector_get(T, idx_E);
		*T_W = contorno(0, y);
	} else if ((i == n - 2) & (j == 0)) {
		*T_N = contorno(x, 0);
		*T_S = gsl_vector_get(T, idx_S);
		*T_E = contorno(1, y);
		*T_W = gsl_vector_get(T, idx_W);
	} else if ((i == 0) & (j == n - 2)) {
		*T_N = gsl_vector_get(T, idx_N);
		*T_S = contorno(x, 1);
		*T_E = gsl_vector_get(T, idx_E);
		*T_W = contorno(0, y);
	} else if ((i == n - 2) & (j == n - 2)) {
		*T_N = gsl_vector_get(T, idx_N);
		*T_S = contorno(x, 1);
		*T_E = contorno(1, y);
		*T_W = gsl_vector_get(T, idx_W);
	} else if (i == 0) {
		*T_N = gsl_vector_get(T, idx_N);
		*T_S = gsl_vector_get(T, idx_S);
		*T_E = gsl_vector_get(T, idx_E);
		*T_W = contorno(0, y);
	} else if (i == n - 2) {
		*T_N = gsl_vector_get(T, idx_N);
		*T_S = gsl_vector_get(T, idx_S);
		*T_E = contorno(1, y);
		*T_W = gsl_vector_get(T, idx_W);
	} else if (j == 0) {
		*T_N = contorno(x, 0);
		*T_S = gsl_vector_get(T, idx_S);
		*T_E = gsl_vector_get(T, idx_E);
		*T_W = gsl_vector_get(T, idx_W);
	} else if (j == n - 2) {
		*T_N = gsl_vector_get(T, idx_N);
		*T_S = contorno(x, 1);
		*T_E = gsl_vector_get(T, idx_E);
		*T_W = gsl_vector_get(T, idx_W);
	} else {
		*T_N = gsl_vector_get(T, idx_N);
		*T_S = gsl_vector_get(T, idx_S);
		*T_E = gsl_vector_get(T, idx_E);
		*T_W = gsl_vector_get(T, idx_W);
	}
}

void calc_jacobiano_analitico(gsl_spmatrix* A, int n, gsl_vector *T) {

	int unknows = (n - 1) * (n - 1);

	int i, j;

	double h = 1.0 / n;

	for (int linha = 0; linha < unknows; linha++) {
		double T_N, T_S, T_E, T_W, T_O;
		define_vizinhos(linha, T, n, &T_O, &T_N, &T_S, &T_E, &T_W);

		i = linha % (n - 1);
		j = linha / (n - 1);
		int linha_N = i + (j - 1) * (n - 1);
		int linha_S = i + (j + 1) * (n - 1);
		int linha_E = i + 1 + j * (n - 1);
		int linha_W = i - 1 + j * (n - 1);

		double k_n = k_interface(T_N, T_O);
		double k_s = k_interface(T_S, T_O);
		double k_e = k_interface(T_E, T_O);
		double k_w = k_interface(T_W, T_O);

		double dk_n = dk_interface(T_N, T_O);
		double dk_s = dk_interface(T_S, T_O);
		double dk_e = dk_interface(T_E, T_O);
		double dk_w = dk_interface(T_W, T_O);

		double deriv = (T_N * dk_n + T_S * dk_s + T_E * dk_e + T_W * dk_w
				- (dk_n + dk_s + dk_e + dk_w) * T_O - (k_n + k_s + k_e + k_w))
				/ (h * h);

		gsl_spmatrix_set(A, linha, linha, deriv);
		if (i > 0)
			gsl_spmatrix_set(A, linha, linha_W,
					(k_w + dk_w * (T_W - T_O)) / (h * h));
		if (i < n - 2)
			gsl_spmatrix_set(A, linha, linha_E,
					(k_e + dk_e * (T_E - T_O)) / (h * h));
		if (j > 0)
			gsl_spmatrix_set(A, linha, linha_N,
					(k_n + dk_n * (T_N - T_O)) / (h * h));
		if (j < n - 2)
			gsl_spmatrix_set(A, linha, linha_S,
					(k_s + dk_s * (T_S - T_O)) / (h * h));

	}
}

void solver_linear( gsl_splinalg_itersolve *solver, gsl_spmatrix* A , gsl_vector* b, gsl_vector* x, float tol, int *iter, double *rnorm){
	size_t max_iter = 2000;
	//Resolucao de Sistema Linear
	*iter = 0;
	int status;

	do{
		status = gsl_splinalg_itersolve_iterate(A, b, tol, x, solver);
		*iter += 1;
	}while(status == GSL_CONTINUE && *iter < max_iter);
	*rnorm = gsl_splinalg_itersolve_normr(solver);
}


int main(int argc, char **argv) {

	int n = atoi(argv[1]);
	

	int unknows = (n - 1) * (n - 1);
	int iters, newton_steps = 0;
	double tol_linear = 1.0e-6; /* solution relative tolerance linear solver */
	     /* maximum iterations */
	double tol_newton = 1e-4;
	const gsl_splinalg_itersolve_type *Solver = gsl_splinalg_itersolve_gmres;
	double rnorm;
	int use_jacobi = 0;
	gsl_splinalg_itersolve *gmres_solver = gsl_splinalg_itersolve_alloc(Solver, unknows, 0);

	if(argc > 2){
		use_jacobi = atoi(argv[2]);
	}

	if(argc > 3){
		tol_newton = atof(argv[3]);
		tol_linear = atof(argv[4]);
	}

	gsl_vector* T = gsl_vector_alloc(unknows);
	for (int i = 0; i < unknows; ++i) {
		gsl_vector_set(T, i, 0);
	}

	gsl_spmatrix* A = gsl_spmatrix_alloc_nzmax(unknows, unknows, 5 * unknows, GSL_SPMATRIX_TRIPLET);
	gsl_spmatrix* C;
	//Alocando Vetores
	gsl_vector* r = gsl_vector_alloc(unknows);
	gsl_vector* dx = gsl_vector_alloc(unknows);
	gsl_vector* jacobi = gsl_vector_alloc(unknows);
	gsl_vector_set_zero (dx);
	printf("*************************************************************************\n");
	printf("n = %d \n", n);
	printf("Toleracia Newton = %e\n", tol_newton);
	printf("Toleracia Linear = %e\n", tol_linear);
	printf("Precondicionador = %d\n", use_jacobi);
	printf("Qtd de Variaveis = %d\n", unknows);

	//Calculo do Residuo Inicial
	calc_residuo(n, T, r);
	gsl_blas_dscal (-1.0, r);
	rnorm = gsl_blas_dnrm2(r);
	printf("Residuo Inicial = %lf\n\n", rnorm);
	printf("*************************************************************************\n");

	int total_iters = 0;

	//Metodo de Newton
	do {
		printf("Passo de Newton %d\n", newton_steps++);

		printf("\tCalculo Jacobiano\n");
		calc_jacobiano(A, n, T);
		printf("\tFim Calculo Jacobiano\n");

		//Aplicar Precondicionador Jacobi
		if(use_jacobi){
			get_jacobi(A, jacobi);
			apply_jacobi_mat(A, jacobi);
			apply_jacobi_vec(r, jacobi);
		}

		C = gsl_spmatrix_crs(A);
		printf("\tInicio Solver Linear\n");
		solver_linear(gmres_solver, C, r, dx, tol_linear, &iters, &rnorm);
		printf("\tFim Solver Linear\n");
		printf("\tIteracoes Lineares = %d\n\trnorm_linear = %e\n", iters, rnorm);
		total_iters += iters;

		//Atualizacao de T
		gsl_blas_daxpy(1.0, dx, T);

		calc_residuo(n, T, r);
		gsl_blas_dscal (-1.0, r);
		rnorm = gsl_blas_dnrm2(r);

		printf("\trnorm_newton = %e\n\n", rnorm);
	} while(rnorm>tol_newton);

	printf("Fim do Método de Newton\n");

	//Escrita da solucao e do arquivo de saida 
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

	double difsol = gsl_blas_dnrm2(sol)/unknows;

	fprintf(saida, "\n\tNorma Residuo = %e\n", difsol);
	printf("\nNorma Residuo = %e\n", difsol);

	fclose(saida);

	printf("%d\t%d\t%d\t%e\n", n, newton_steps, total_iters, difsol/unknows);

	gsl_vector_free(sol);
	gsl_vector_free(T);
	gsl_vector_free(r);	
	gsl_vector_free(dx);	
	gsl_vector_free(jacobi);
	gsl_splinalg_itersolve_free(gmres_solver);

}
