#include <stdio.h>
#include <math.h>
#include "linearsys.h"

void qr_test( char * file_name, int nrows, int ncols )
{
	printf("Testing QR decomposition of matrix A:\n");
	Matrix A = mat_new_from_file( file_name, nrows, ncols );

	mat_print_stdout(&A, "%g", true);

	Matrix Q = {0};
	Matrix R = {0};
	linst_qr_mgs(&A, &Q, &R);
	
	printf("We found\nQ:\n");
	mat_print_stdout(&Q,"%g", true);
	printf("R:\n");
	mat_print_stdout(&R,"%g", true);

	Matrix QR = mat_mult_new(&Q,&R);
	Matrix QONC = {0};
	mat_gram(&Q, &QONC);
	
	printf("\nProduct QR is:\n");
	mat_print_stdout(&QR,"%g", true);
	
	printf("\nProduct QtQ is:\n");
	mat_print_stdout(&QONC,"%g", true);
	
	printf("----------------------------------------\n\n");
	
	Matrix *to_free[]={&Q,&R,&A, &QR, &QONC};
	mat_free_many(to_free, 5);
	
}

int main()
{
	
	qr_test( "samplematrix1.txt", 3, 2 );
	qr_test( "samplematrix2.txt", 3, 3 );
	qr_test( "samplematrix3.txt", 4, 3 );
	
	return 0;
}

