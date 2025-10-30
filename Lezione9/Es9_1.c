#include "comp_physics.h"
#include <stdio.h>
#include <math.h>

void print_qr_mats( MatrixDouble A, QR_Mats qr, MatrixDouble QR, MatrixDouble QONC)
{
	printf("Starting matrix:\n");
	printMatDGraph( A );
	printf("QR product is:\n");
	printMatDGraph(QR);
	printf("Matrix Q is:\n");
	printMatDGraph( qr.Q );
	printf("Matrix R is:\n");
	printMatDGraph( qr.R );
	printf("QTQ is:\n");
	printMatDGraph( QONC );
	putchar('\n');
}
// possiamo calcolare la bontà del risultato confrontando la differenza fra
// A e QR e il condition number

MatrixDouble test_ONC( MatrixDouble * Q )
{
	MatrixDouble tQ = allocMatD( Q->ncols, Q->nrows );
	transposeMatD( *Q, &tQ );
	MatrixDouble id = matMultD( tQ, *Q, NULL );
	freeMatD( tQ );
	return id;
}

void test_on_A( MatrixDouble A )
{
	QR_Mats qr = QR_decomp( A );
	MatrixDouble QR = matMultD( qr.Q, qr.R, NULL );
	
	MatrixDouble QONC = test_ONC( &qr.Q );
	
	print_qr_mats( A, qr, QR, QONC );
	free_MatD(qr.Q, qr.R, QR, QONC);
}

void qr_test( char * file_name, int nrows, int ncols )
{
	printf("Testing of %s\n", file_name);
	MatrixDouble A = readMatD( file_name, nrows, ncols );
	test_on_A(A);
}

int main()
{
	MatrixDouble id = allocMatD(4, 4);
	diagMatD( id, 1 );
	test_on_A(id);
	freeMatD(id);
	
	qr_test( "data/samplematrix1.txt", 3, 2 );
	qr_test( "data/samplematrix2.txt", 3, 3 );
	qr_test( "data/samplematrix3.txt", 4, 3 );
	
	eprint("program ended fine");
	
	return 0;
}

