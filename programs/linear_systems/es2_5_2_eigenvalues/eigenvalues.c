#include "linearsys.h"
#include <stdio.h>

void test_pure_qr(const char *filename, int nrows, int ncols, double tol ){
	
	Matrix test = mat_new_from_file(filename,nrows,ncols);
	Array eigenvalues = {0};
	Matrix eigenvectors = {0};
	int num_of_iterations=0;
	
	printf("Matrix A to decompose:\n");
	mat_print_stdout(&test, "%g",true);
	
	num_of_iterations=linst_pure_qr(&test, &eigenvalues, &eigenvectors, tol);
	
	printf("Eigenvalues:\n");
	arr_print_inline(&eigenvalues,"% g",true, true);
	printf("Eigenvector matrix:\n");
	mat_print_stdout(&eigenvectors, "% 3.4f", true);
	
	printf("For each pair of eigenvectors X and eigenvalues λ found,"
		"the following are the differences between λX and AX:\n");
	
	Array evec={0};
	Array evproduct={0};
	for(int i=0; i<ncols; i++)
	{
		mat_col_to_arr(&eigenvectors, &evec,i);
		mat_vecmult(&test, &evec, &evproduct);
		arr_axpy(&evproduct, -ARR(eigenvalues,i), &evec);
		arr_print_inline(&evproduct,"%g",true,true);
	}
	
	printf("Number of iterations: %d\n", num_of_iterations);
	
	arr_free(&eigenvalues);
	mat_free(&eigenvectors);
	mat_free(&test);
	
}

int main()
{
	test_pure_qr("samplematrix.txt",4,4, 1e-10);
	
	exit(EXIT_SUCCESS);
}
