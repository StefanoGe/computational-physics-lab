#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "linearsys.h"

Matrix test_matrix( double eps, bool swap_rows)
{
	Matrix A = mat_new(2,2);
	MAT(A,0,0) = -eps;
	MAT(A,1,0) = 1;
	MAT(A,0,1) = 1;
	MAT(A,1,1) = -1;
	
	if(swap_rows)
		mat_swap_row(&A,0,1);
	
	return A;
}

void cond_numb_study( double eps, bool swap_rows )
{
	printf("Testing for eps = %g\n", eps);
	Matrix A = test_matrix(eps, swap_rows);
	
	double x_[]={1,1};
	Array x = arr_asarr(x_,2);
	Array b = mat_vecmult_new(&A,&x);
	
	printf("A is: \n");
	mat_print_stdout(&A,"%g",true);
	
	printf("x is: ");
	arr_print_inline( &x, "%g", true, true );
	
	printf("b=Ax is: ");
	arr_print_inline( &b, "%.30g", true, true );

	Matrix LU={0};
	mat_cp(&A, &LU);
	linst_lu_factor_no_pivot(&LU, 0);
	Matrix L=linst_lu_extract_l(&LU);
	Matrix U=linst_lu_extract_u(&LU);

	printf("L is:\n");
	mat_print_stdout( &L, "%g", true );
	
	printf("U is:\n");
	mat_print_stdout( &U, "%g", true );

	linst_forwsubst_inplace( &LU, &b, true );
	linst_backsubst_inplace( &LU, &b );
	
	printf("x found by solving LUx=b: \n");
	arr_print_inline(&b, "%g", true, true);
	
	printf("x_calc - x_true is: ");
	arr_axpy(&b, -1, &x);
	arr_print_inline( &b, "%g", true, true );
	
	Matrix LU_product = mat_mult_new( &L, &U );
	
	mat_axpy(&A, -1, &LU_product);
	
	printf( "A - LU:\n");
	mat_print_stdout( &A, "%g", true );

	
	
	
	printf("------------------------------\n");
}



int main()
{
	//test_inv();
	cond_numb_study( 1e-12, false );
	cond_numb_study( 1e-20, false );
	cond_numb_study( 1e-20, true );
	
	
	return 0;
}

