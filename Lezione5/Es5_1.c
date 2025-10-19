// Es5_1.c LU Decomposition

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "DoubleUtilities.c"
#include "subst"

#define PI 3.1415926535897932384


void testLU( char * filename, int dim )
{
	MatrixDouble testMat = readMatD( filename, dim, dim);
	LUMats lumats = LUDecompLow ( testMat );
	
	printf("Test of %s:\n", filename);
	printMatDGraph(testMat);
	printf("Upper is: \n");
	printMatDGraph( lumats.U );
	printf("Lower is: \n");
	printMatDGraph( lumats.L );
	
	freeMatD(testMat);
	freeMatD(lumats.U);
	freeMatD(lumats.L);
}

MatrixDouble planeTranslation( double x, double y )
{
	MatrixDouble matT = allocMatD( 3, 3 );
	diagMatD( matT, 1 );
	
	matT.val[0][2] = x;
	matT.val[1][2] = y;
	
	printMatDGraph(matT);
	
	return matT;
}

MatrixDouble planeRotation( double theta )
{
	MatrixDouble matR = allocMatD( 3, 3 );
	
	matR.val[0][0] = cos(theta);
	matR.val[1][1] = matR.val[0][0];
	matR.val[2][2] = 1;
	matR.val[0][1] = sin(theta);
	matR.val[1][0] = - matR.val[0][1];
	
	printMatDGraph(matR);\
	
	return matR;
}

void testRotTra()
{
	MatrixDouble Tleft = planeTranslation(3, -1);
	MatrixDouble R = planeRotation( PI/5 );
	MatrixDouble Tright = planeTranslation(-3, 1);
	
	ArrayDouble z = buildArrD( 3, 2.0, 2.0, 1.0 );
	
	printf("z is\n");
	printArrDPar(z);
	
	MatrixDouble A = matMultD( Tleft, R, CREATE_MAT );
	matMultD( A, Tright, &A );
	
	printf("\nA matrix is:\n");
	printMatDGraph( A );
	
	MatrixDouble z_mat = asMatrix( z );
	MatrixDouble b = matMultD( A, z_mat, CREATE_MAT );
	
	printf("b arr is: \n");
	printMatDGraph( b );
	
	LUMats lumat = LUDecompLow( A );
	printf("Upper is: \n");
	printMatDGraph( lumat.U );
	printf("Lower is: \n");
	printMatDGraph( lumat.L );
	
	ArrayDouble b_arr = allocArrD( b.nrows );
	for(int i = 0; i < b.nrows; i++)
		b_arr.val[i] = b.val[i][0];
	
	ArrayDouble y = forwSubst( lumat.L, b_arr );
	ArrayDouble x = backSubst( lumat.U, y );
	
	printf("y is: \n");
	printArrDPar( y );
	printf("x is: \n");
	printArrDPar( x );
	
	freeArrD(z);
	freeArrD(y);
	freeArrD(x);
	
	freeAllMatD( R, Tleft, Tright, z_mat, b, A, NULL_MAT );
}

int main()
{
	testLU( "samplematrix1.txt", 3 );
	testLU( "samplematrix2.txt", 4 );
	testLU( "samplematrix3.txt", 4 );
	
	testRotTra();
	
	
	return 0;
}

/*
void testDiag()
{
	MatrixDouble mat = allocMatD(5, 5);
	diagMatD( mat, -3.7 );
	printMatDGraph(mat);
	freeMatD( mat );
}

void testSlice()
{
	MatrixDouble testMat = readMatD("samplematrix1.txt", 3, 3);
	MatrixDouble slices = sliceMatD( testMat, "1 2", "1" );
	printf("Sliced matrix is:\n");
	printMatDGraph(slices);
	
	freeMatD( testMat );
	freeMatD( slices );
}
*/
