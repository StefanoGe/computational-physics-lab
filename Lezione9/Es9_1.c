#include "comp_physics.h"
#include <stdio.h>
#include <math.h>

typedef struct{
	MatrixDouble Q;
	MatrixDouble R;
} QR_Mats;

QR_Mats alloc_QR( int nrows, int ncols  )
{
	QR_Mats qr_mats;
	qr_mats.Q = allocMatD( nrows, ncols );
	qr_mats.R = allocMatD( ncols, ncols );
	return qr_mats;
}

double norm( MatrixDouble * Q, int col )
{
	double norm_sq = 0;
	for( int i = 0; i < Q->nrows; i++ )
		norm_sq += Q->val[i][col] * Q->val[i][col];
		
	return sqrt(norm_sq);
}

double normalise( MatrixDouble * Q, int col )
{		
	double vec_norm = norm( Q, col );
	for( int i = 0; i < Q->nrows; i++ )
		Q->val[i][col] /= vec_norm;
		
	return vec_norm;
}

QR_Mats QR_decomp( MatrixDouble A )
{
	const int height = A.nrows;
	const int nvecs = A.ncols;
	QR_Mats qr = alloc_QR( height, nvecs );
	diagMatD(qr.R, 0);
	cpMatToMatD( A, qr.Q, 0, 0 );
	
	// va aggiunta una cosa preliminare per avere q0 (gli altri q hanno bisogno 
	// di lui)
	
	// loop che controlla il q che sto usando
	for( int i = 0; i < nvecs; i ++ )
	{
		qr.R.val[i][ i ] = normalise( &qr.Q, i );
		// Mi serve ora un algoritmo che calcoli i vari vj
		for( int j = i + 1; j < nvecs; j++)
		{
			double scalar_product = 0;
			for( int p = 0; p < height; p++ )
				scalar_product += qr.Q.val[p][i] * qr.Q.val[p][j];
			
			qr.R.val[i][j] = scalar_product;
			
			for( int p = 0; p < height; p++ )
			{
				qr.Q.val[p][j] -= scalar_product * qr.Q.val[p][i];
			}
		}
	}
	// freeAllMatD( NULL_MAT );
	return qr;
}

void print_qr_mats( MatrixDouble A, QR_Mats qr)
{
	printf("Starting matrix:\n");
	printMatDGraph( A );
	printf("Matrix Q is:\n");
	printMatDGraph( qr.Q );
	printf("Matrix R is:\n");
	printMatDGraph( qr.R );
	putchar('\n');
}

void test_on_A( MatrixDouble A )
{
	QR_Mats qr = QR_decomp( A );
	print_qr_mats( A, qr );
	freeAllMatD(qr.Q, qr.R, NULL_MAT);
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

