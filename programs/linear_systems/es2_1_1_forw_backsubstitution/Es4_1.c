#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "genutil.h"
#include "linearsys.h"

void testID()
{
	Matrix identity = mat_new(3,3);
	mat_diag( &identity, 1 );
	double b[] = {3.14, -1.0, 2.1};
	Array b_arr = arr_asarr(b, 3);
	Array sol = {0};
	linst_forwsubst(&identity, &b_arr, &sol, false);
	Array sol_back = {0};
	linst_backsubst(&identity, &b_arr, &sol_back);
	
	printf("Test of forward and back substitution on the identity matrix.\n"
			"Starting right-hand side of Ix=b and expected solution:\n");
	arr_print_inline(&sol, "%g", true, true);
	printf("Forward substitution\n");
	arr_print_inline( &sol, "%g", true, true );
	printf("Back-substitution:\n");
	arr_print_inline( &sol_back, "%g", true, true);
	
	mat_free(&identity);
	arr_free(&sol);
	arr_free(&sol_back);
}

void test_forward()
{
	Matrix test_mat1 = mat_new_from_file( "samplematrix1.txt", 3, 3 );
	double b1[] = {-4.0, 2.0, 1.0};
	Array b1_arr = arr_asarr(b1,3);
	Array sol1 = {0};
	linst_forwsubst( &test_mat1, &b1_arr, &sol1, false );
	
	Matrix test_mat2 = mat_new_from_file( "samplematrix2.txt", 4, 4 );
	double b2[] = {-4.0, 1.0, -3.0, 5.0};
	Array b2_arr = arr_asarr(b2, 4);
	Array sol2 = {0};
	linst_forwsubst( &test_mat2, &b2_arr, &sol2, false );
	
	// Print of 1
	
	puts("Test of forward substitution in two special cases, to solve\n"
				"the linear system Lx=b\n");
	
	puts("First matrix L:");
	
	mat_print_stdout(&test_mat1, "%2g", true);
	
	puts("Right-hand side b:");
	
	carr_print_inline(b1, 3, "%g", true, true);
	
	puts("\nSolution X:");
	
	arr_print_inline(&sol1, "%g", true, true);
	
	// Print of 2
	
	printf("\nSecond matrix L:\n");
	
	mat_print_stdout(&test_mat2, "%2g", true);
	
	puts("Right-hand side b:");
	
	carr_print_inline(b2, 4, "%g", true, true);
	
	puts("\nSolution X:");
	
	arr_print_inline(&sol2, "%g", true, true);
	
	Matrix *mats[] = {&test_mat1, &test_mat2};
	mat_free_many(mats,2);
	Array *arrs[]={&sol1, &sol2};
	arr_free_many(arrs, 2);
}

void test_back()
{
	Matrix test_mat1 = mat_new_from_file( "samplematrix1back.txt", 3, 3 );
	double b1[] = {1, 1, 6};
	Array b1_arr = arr_asarr(b1, 3);
	Array sol1 = {0};
	linst_backsubst( &test_mat1, &b1_arr, &sol1 );
	
	Matrix test_mat2 = mat_new_from_file( "samplematrix2back.txt", 4, 4 );
	double b2[] = {4.0, 1.0, 1.0, 5.0};
	Array b2_arr = arr_asarr(b2, 4);
	Array sol2 = {0};
	linst_backsubst( &test_mat2, &b2_arr, &sol2 );
	
	// Print of 1
	
	puts("Test of backward substitution in two special cases, to solve\n"
				"the linear system UX=b\n");
	
	puts("First matrix L:");
	
	mat_print_stdout(&test_mat1, "%2g", true);
	
	puts("Right-hand side b:");
	
	carr_print_inline(b1, 3, "%g", true, true);
	
	puts("\nSolution X:");
	
	arr_print_inline(&sol1, "%g", true, true);
	
	// Print of 2
	
	printf("\nSecond matrix U:\n");
	
	mat_print_stdout(&test_mat2, "%2g", true);
	
	puts("Right-hand side b:");
	
	carr_print_inline(b2, 4, "%g", true, true);
	
	puts("\nSolution X:");
	
	arr_print_inline(&sol2, "%g", true, true);
	
	Matrix *mats[] = {&test_mat1, &test_mat2};
	mat_free_many(mats,2);
	Array *arrs[]={&sol1, &sol2};
	arr_free_many(arrs, 2);
}

int main()
{
	testID();
	puts("\n-------------------------------\n");
	test_forward();
	puts("\n-------------------------------\n");
	test_back();
	
	return 0;
}
