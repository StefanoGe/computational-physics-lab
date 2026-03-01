
#ifdef HAVE_LAPACKE

#include <lapacke.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "genutil.h"
#include "linearsys.h"

void test_forward()
{
	
	Matrix test_mat1 = mat_new_from_file( "samplematrix1.txt", 3, 3 );
	double b1[] = {-4.0, 2.0, 1.0}; // gets overwritten by lapack
	
	// Print of 1
	
	puts("Test of forward substitution in two special cases, to solve\n"
				"the linear system Lx=b with LAPACK\n");
	
	puts("First matrix L:");
	
	mat_print_stdout(&test_mat1, "%2g", true);
	
	puts("Right-hand side b:");
	
	carr_print_inline(b1, 3, "%g", true, true);
	
	puts("\nSolution X:");
	
	int info1 = LAPACKE_dtrtrs(LAPACK_ROW_MAJOR, 
		'L', 
		'N',
		'N', 
		3, 
		1, 
		test_mat1.data, 
		test_mat1.ncols, 
		b1, 
		1 );
	
	if(info1 != 0)
		printf("lapack info: %d\n\n", info1);
	
	carr_print_inline(b1, 3, "%g", true, true);
	
	Matrix test_mat2 = mat_new_from_file( "samplematrix2.txt", 4, 4 );
	double b2[] = {-4.0, 1.0, -3.0, 5.0};
	
	// Print of 2
	
	printf("\nSecond matrix L:\n");
	
	mat_print_stdout(&test_mat2, "%2g", true);
	
	puts("Right-hand side b:");
	
	carr_print_inline(b2, 4, "%g", true, true);
	
	puts("\nSolution X:");
	
	int info2 = LAPACKE_dtrtrs(LAPACK_ROW_MAJOR, 
		'L', 
		'N',
		'N', 
		test_mat2.nrows, 
		1, 
		test_mat2.data, 
		test_mat2.ncols, 
		b2, 
		1 );
	
	if(info2 != 0)
		printf("lapack info: %d\n\n", info2);
	
	carr_print_inline(b2, 4, "%g", true, true);
	
}

void test_back()
{
	Matrix test_mat1 = mat_new_from_file( "samplematrix1back.txt", 3, 3 );
	double b1[] = {1, 1, 6};
	
	// Print of 1
	
	puts("Test of backward substitution in two special cases, to solve\n"
				"the linear system UX=b\n");
	
	puts("First matrix L:");
	
	mat_print_stdout(&test_mat1, "%2g", true);
	
	puts("Right-hand side b:");
	
	carr_print_inline(b1, 3, "%g", true, true);
	
	puts("\nSolution X:");
	
	int info1 = LAPACKE_dtrtrs(LAPACK_ROW_MAJOR, 
		'U', 
		'N',
		'N', 
		3, 
		1, 
		test_mat1.data, 
		test_mat1.ncols, 
		b1, 
		1 );
	
	if(info1 != 0)
		printf("lapack info: %d\n\n", info1);
	
	carr_print_inline(b1, 3, "%g", true, true);
	
	Matrix test_mat2 = mat_new_from_file( "samplematrix2back.txt", 4, 4 );
	double b2[] = {4.0, 1.0, 1.0, 5.0};
	Array sol2 = {0};
	linst_backsubst( test_mat2.rows, 4, b2, &sol2 );
	
	// Print of 2
	
	printf("\nSecond matrix U:\n");
	
	mat_print_stdout(&test_mat2, "%2g", true);
	
	puts("Right-hand side b:");
	
	carr_print_inline(b2, 4, "%g", true, true);
	
	puts("\nSolution X:");
	
	int info2 = LAPACKE_dtrtrs(LAPACK_ROW_MAJOR, 
		'U', 
		'N',
		'N', 
		test_mat2.nrows, 
		1, 
		test_mat2.data, 
		test_mat2.ncols, 
		b2, 
		1 );
	
	if(info2 != 0)
		printf("lapack info: %d\n\n", info2);
	
	carr_print_inline(b2, 4, "%g", true, true);
}

int main()
{
	puts("\n-------------------------------\n");
	test_forward();
	puts("\n-------------------------------\n");
	test_back();
	
	return 0;
}



#endif


#ifndef HAVE_LAPACKE
int main(){return 1;}
#endif

/*
// check_lapacke.c
// Compile example: gcc -o check_lapacke check_lapacke.c -llapacke -llapack -lblas -lgfortran -lm

#include <stdio.h>

#ifdef HAVE_LAPACKE
#  include <lapacke.h>
#else
#  error "HAVE_LAPACKE is not defined — add -DHAVE_LAPACKE=1 to compiler flags"
#endif

int main(void) {
    printf("LAPACKE header is included.\n");

    // Very minimal runtime test: check one symbol from lapacke.h
    // We call a function that exists in almost every LAPACKE installation
    int major, minor, patch;
    LAPACKE_ilaver(&major, &minor, &patch);

    printf("LAPACKE is linked and usable.\n");
    printf("LAPACK version: %d.%d.%d\n", major, minor, patch);

    // Optional: test one solver routine (dgetrf – LU decomposition)
    const int n = 2;
    double A[4] = {4.0, 3.0, 6.0, 3.0};   // row-major
    int ipiv[2];
    int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv);

    if (info == 0) {
        printf("LAPACKE_dgetrf works (simple 2×2 test passed).\n");
    } else {
        printf("LAPACKE_dgetrf failed (info = %d)\n", info);
    }

    return 0;
}
*/
