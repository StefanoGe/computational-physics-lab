// Es6_1.c LU Decomposition

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearsys.h"
#include "genutil.h"


void testLUP( char * filename, int dim )
{
	Matrix test_mat = mat_new_from_file(filename, dim, dim);
	int *pivots = malloc(dim*sizeof(int));
	int info_fail;
	
	printf("Testing the lu decomposition of matrix %s:\n", filename);
	mat_print_stdout(&test_mat, "%3g ", true);
	
	if( ( info_fail=linst_lup_factor(&test_mat, pivots, 0) ) )
	{
		eprint("factorization failed at n = %d",info_fail);
		return;
	}
	
	Matrix U = linst_lu_extract_u(&test_mat);
	Matrix L = linst_lu_extract_l(&test_mat);
	Matrix P = linst_lu_extract_p(pivots, dim);
	
	printf("Upper is: \n");
	mat_print_stdout(&U, "%3g ", true);
	
	printf("Lower is: \n");
	mat_print_stdout(&L, "%3g ", true);
	
	printf("P is:\n");
	mat_print_stdout(&P, "%3g ", true);
	
	printf("Determinant is: %g\n", linst_lup_det(&test_mat, pivots));
	printf("---------------------------------------\n");
	
	Matrix *to_free[]={&test_mat, &U, &L, &P};
	mat_free_many(to_free, 4);
	free(pivots);
}

Matrix planeTranslation( double x, double y )
{
	Matrix T = mat_new( 3, 3 );
	mat_diag(&T, 1);
	
	MAT(T,0,2) = x;
	MAT(T,1,2) = y;
	
	return T;
}

Matrix planeRotation( double theta )
{
	Matrix R = mat_new( 3, 3 );
	
	MAT(R,0,0) = cos(theta);
	MAT(R,1,1) = MAT(R,0,0);
	MAT(R,2,2) = 1;
	MAT(R,0,1) = sin(theta);
	MAT(R,1,0) = - MAT(R,0,1);
	
	return R;
}

Matrix TR_composition(  )
{
	Matrix Tleft = planeTranslation(3, -1);
	Matrix R = planeRotation( PI/5 );
	Matrix Tright = planeTranslation(-3, 1);
	
	Matrix A = mat_mult_new( &Tleft, &R );
	Matrix B = mat_mult_new( &A, &Tright );
	
	Matrix *to_free[]={&Tleft, &R, &Tright, &A};
	mat_free_many(to_free, 4);
	
	return B;
}

void testRotTra()
{
	Matrix A = TR_composition();
	int *pivots=malloc(sizeof(int)*3);
	double z[]={2.0, 2.0, 1.0};
	Array z_arr=arr_asarr(z,3);
	
	printf("\nA matrix is:\n");
	mat_print_stdout( &A, "%g", true );
	
	printf("z is\n");
	arr_print_inline(&z_arr,"%g",true,true);
	
	Array b = mat_vecmult_new( &A, &z_arr );
	
	printf("b = Az is: \n");
	arr_print_inline(&b,"%g",true,true);
	
	linst_lup_factor(&A,pivots,0);
	Matrix U=linst_lu_extract_u(&A);	
	Matrix L=linst_lu_extract_l(&A);
	Matrix P=linst_lu_extract_p(pivots,3);

	printf("P is: \n");
	mat_print_stdout( &P, "%g", true );
	printf("LU decomposition\nU is: \n");
	mat_print_stdout( &U, "%g", true );
	printf("L is: \n");
	mat_print_stdout( &L, "%g", true );
	
	linst_lup_solve_inplace(&A, pivots,&b);
	
	printf("x of Ax = b is: \n");
	arr_print_inline( &b, "%g",true,true );
	
	arr_free(&b);
	arr_free(&z_arr);
	free(pivots);
	Matrix *mat_to_free[]={&U,&L,&A};
	mat_free_many(mat_to_free,3);
}


int main()
{

	testLUP( "samplematrix1.txt", 3 );
	testLUP( "samplematrix2.txt", 4 );
	testLUP( "samplematrix3.txt", 4 );

	testRotTra();
	
	return 0;
}
