#include "linearsys.h"
#include <stdio.h>

void test_pure_qr(const char *filename, int nrows, int ncols, double tol, 
	int max_iter){
	
	Matrix test = mat_new_from_file(filename,nrows,ncols);
	Array eigenvalues = {0};
	Matrix eigenvectors = {0};
	
	
	printf(ANSI_YELLOW"Matrix A to decompose:\n"ANSI_RESET);
	mat_print_stdout(&test, "%g",true);
	
	if(linst_pure_qr(&test, &eigenvalues, &eigenvectors, 
		tol, max_iter))
		raiseErr("Pure QR did not reach precision of %g after %d\n"
				 "iterations", tol, max_iter);
	
	printf(ANSI_YELLOW"Eigenvalues:\n"ANSI_RESET);
	arr_print_inline(&eigenvalues,"% g",true, true);
	printf(ANSI_YELLOW"Eigenvector matrix:\n"ANSI_RESET);
	mat_print_stdout(&eigenvectors, "% 3.4f", true);
	
	printf(ANSI_YELLOW"For each pair of eigenvectors X and eigenvalues λ found,\n"
		"the following are the differences between λX and AX:\n"ANSI_RESET);
	
	Array evec={0};
	Array evproduct={0};
	for(int i=0; i<ncols; i++)
	{
		printf("X%d: ",i+1);
		mat_col_to_arr(&eigenvectors, &evec,i);
		mat_vecmult(&test, &evec, &evproduct);
		arr_axpy(&evproduct, -ARR(eigenvalues,i), &evec);
		arr_print_inline(&evproduct,"%+8.2g",true,true);
	}
	
//	printf("Number of iterations: %d\n", num_of_iterations);
	
	arr_free(&eigenvalues);
	mat_free(&eigenvectors);
	mat_free(&test);
	
}

int main()
{
	test_pure_qr("samplematrix.txt",4,4, 1e-10, 1000);
	
	exit(EXIT_SUCCESS);
}
