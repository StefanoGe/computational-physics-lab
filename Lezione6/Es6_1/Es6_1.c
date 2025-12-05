// Es6_1.c LU Decomposition

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comp_physics.h"

#define PI 3.1415926535897932384


void testLUP( char * filename, int dim )
{
	MatrixDouble testMat = readMatD( filename, dim, dim);
	
	LUP_Mats lup_mats = LUP_decomposition ( testMat );
	
	printf("Test of %s:\n", filename);
	printMatDGraph(testMat);
	
	printf("Upper is: \n");
	printMatDGraph( lup_mats.U );
	
	printf("Lower is: \n");
	printMatDGraph( lup_mats.L );
	
	printf("P is:\n");
	printMatDGraph(lup_mats.P);
	
	printf("Determinant is: %lf\n", detTri(lup_mats.U));
	printf("---------------------------------------\n");
	
	freeAllMatD(testMat, lup_mats.U, lup_mats.P, lup_mats.L, NULL_MAT);
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

MatrixDouble buildAMatrix(  )
{
	MatrixDouble Tleft = planeTranslation(3, -1);
	MatrixDouble R = planeRotation( PI/5 );
	MatrixDouble Tright = planeTranslation(-3, 1);
	
	MatrixDouble A = matMultD( Tleft, R, CREATE_MAT );
	matMultD( A, Tright, &A );
	
	freeAllMatD( R, Tleft, Tright, NULL_MAT );
	
	return A;
}

void printDebugRotra(MatrixDouble A, MatrixDouble z, MatrixDouble b,
				LUP_Mats lumat, ArrayDouble x, ArrayDouble y)
{
	printf("\nA matrix is:\n");
	printMatDGraph( A );
	
	printf("z is\n");
	printMatDGraph(z);
	
	printf("b arr is: \n");
	printMatDGraph( b );

	printf("Upper is: \n");
	printMatDGraph( lumat.U );
	printf("Lower is: \n");
	printMatDGraph( lumat.L );

	printf("y is: \n");
	printArrDPar( y );
	printf("x is: \n");
	printArrDPar( x );

}

void testRotTra(bool debug)
{
	MatrixDouble A = buildAMatrix();
	MatrixDouble z = asColMatrix(
		buildArrD( 3, 2.0, 2.0, 1.0 ), true);
	
	MatrixDouble b = matMultD( A, z, CREATE_MAT );
	ArrayDouble b_arr = allocArrD( b.nrows );
	
	LUP_Mats lumat = LUP_decomposition( A );
	
	for(int i = 0; i < b.nrows; i++)
		b_arr.val[i] = b.val[i][0];
	
	ArrayDouble y = forwSubst( lumat.L, b_arr );
	ArrayDouble x = backSubst( lumat.U, y );
	
	if(debug) printDebugRotra( A, z, b, lumat, x, y);
	
	freeAllArrD(x, y, b_arr, NULL_ARR);
	
	freeAllMatD( z, b, A, NULL_MAT );
}

int main()
{
	testLUP( "data/samplematrix1.txt", 3 );
	testLUP( "data/samplematrix2.txt", 4 );
	testLUP( "data/samplematrix3.txt", 4 );
	
	testRotTra(true);
	
	
	return 0;
}
