all: linear naolinear naolinearpontofixo transiente

linear: Laplaciano.c solvers.c
	icc -o Laplaciano Laplaciano.c solvers.c -mkl -std=c99

naolinear: LaplacianoNaoLinear_gsl.c
	gcc -o LaplacianoNaoLinear LaplacianoNaoLinear_gsl.c -I/usr/local/include -L/usr/local/lib -lgsl -std=c99

naolinearpontofixo: LaplacianoNaoLinearPontoFixo.c
	icc -o LaplacianoNaoLinearPontoFixo LaplacianoNaoLinearPontoFixo.c solvers.c -mkl -std=c99

transiente: LaplacianoTransiente.c solvers.c
	icc -o LaplacianoTransiente LaplacianoTransiente.c solvers.c -mkl -std=c99
