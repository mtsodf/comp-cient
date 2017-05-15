#ifndef SOME_HEADER_GUARD_WITH_UNIQUE_NAME
#define SOME_HEADER_GUARD_WITH_UNIQUE_NAME

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
#endif

#define TRUE 1
#define FALSE 0
#define MAXITER 2000

typedef struct struct_matriz matriz;

void cg(matriz, double*, double*, double*);
void cg_precon(matriz, double*, double*, double*);
void mat_mult(int, matriz, double* , double*);
int solve_steepest_descent(int, matriz, double*, double*, double, int*);
