#include <stdlib.h>
#include <stdio.h>
#include "comp_physics.h"

MatrixDouble build_matrix( double epsilon)
{
	MatrixDouble mat = allocMatD( 2,2 );
	mat.val[0][0] = -epsilon;
	mat.val[0][1] = 1;
	mat.val[1][0] = 1;
	mat.val[1][1] = -1;
	
	return mat;
}


void cond_numb_study( double eps )
{
	printf("eps = %e\n", eps);
	MatrixDouble mat = build_matrix(eps);
	ArrayDouble x = buildArrD(2, 1.0, 0.0);
	ArrayDouble b = mat_vec_mult( mat, x, NULL );
	
	
	printf("A is: ");
	printMatDGraph(mat);
	printf("x is: ");
	printArrDPar( x, "%.30lf " );
	printf("b is: ");
	printArrDPar( b, "%.30lf " );

	LUMats lumat= LUDecompLow( mat );
	ArrayDouble z = forwSubst( lumat.L, b );
	ArrayDouble x_calc = backSubst( lumat.U, z );
	
	printArrDPar( x_calc, "%.30lf " );
	
	MatrixDouble LU_product = matMultD( lumat.L, lumat.U, CREATE_MAT );
	
	MatrixDouble diff_matrix = mat_diffD( mat, LU_product, CREATE_MAT);
	
	printf( "A - LU:\n");
	printMatDGraph( diff_matrix );
	
	// Condition number in sup-norm
	
	
	
	printf("------------------------------\n");
}

void test_inv()
{
	MatrixDouble id = readMatD( "data/inv_test_matrix.txt", 3, 3 );
	MatrixDouble known_terms = readMatD( "data/inv_test_matrix_2.txt", 3, 3 );
	
	MatrixDouble inv_id = mat_inv(id);
	
	eprint("ok 4\n");
	
	MatrixDouble inv_other = mat_inv(known_terms);
	
	printf("inv id:\n");
	printMatDGraph(inv_id);
	printf("Other inv:\n");
	printMatDGraph(inv_other);
}

/*
	 cond_numb_study( 1e-12 );
	 cond_numb_study( 1e-20 );
*/

int main()
{
	test_inv();

	
	return 0;
}

