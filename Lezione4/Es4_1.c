#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "DoubleUtilities.c"

void checkMatrixColHeight(MatrixDouble matrix, ArrayDouble colVec)
{
	if( matrix.nrows != colVec.length )
	{
		printf("Matrix and colVec should have same height! err in func forwSubst.\n");
		exit(EXIT_FAILURE);
	}
	if( matrix.ncols != matrix.nrows )
	{
		printf("La matrice deve essere quadrata! err in func forwSubst.\n");
		exit(EXIT_FAILURE);
	}
}

ArrayDouble forwSubst( MatrixDouble matrix, ArrayDouble colVec)
{
	checkMatrixColHeight( matrix, colVec );
	
	const int n = colVec.length;
	
	ArrayDouble solutions = allocArrD( n );
	
	double partSum;
	
	for( int row = 0; row < n; row++ )
	{
		partSum = 0;
		for( int col = 0; col < row; col++ )
			partSum += solutions.val[col] * matrix.val[row][col];
		solutions.val[row] = (colVec.val[row] - partSum)/matrix.val[row][row];
	}
	
	return solutions;
}

ArrayDouble backSubst( MatrixDouble matrix, ArrayDouble colVec)
{
	checkMatrixColHeight( matrix, colVec );
	
	const int n = colVec.length;
	
	ArrayDouble solutions = allocArrD( n );
	
	double partSum;
	
	for( int row = n - 1; row >= 0; row-- )
	{
		partSum = 0;
		for( int col = row + 1; col < n; col++ )
			partSum += solutions.val[col] * matrix.val[row][col];
		solutions.val[row] = (colVec.val[row] - partSum)/matrix.val[row][row];
	}
	
	return solutions;
}

void testID()
{
	MatrixDouble testMatrix = readMatD( "startmatrix.txt", 3, 3 );
	ArrayDouble colVec = buildArrD(3, 3.14, -1.0, 2.1);
	ArrayDouble testSol = forwSubst( testMatrix, colVec );
	ArrayDouble testSolBack = backSubst( testMatrix, colVec );
	
	printArrDPar( testSol );
	printArrDPar( testSolBack );
	
	freeArrD(colVec);
	freeMatD( testMatrix );
	freeArrD( testSol );
	freeArrD( testSolBack );
}

void testSample1()
{
	MatrixDouble testMatrix = readMatD( "samplematrix1.txt", 3, 3 );
	ArrayDouble colVec = buildArrD(3, -4.0, 2.0, 1.0);
	ArrayDouble testSol = forwSubst( testMatrix, colVec );
	
	printArrDPar( testSol );
	
	freeArrD(colVec);
	freeMatD( testMatrix );
	freeArrD( testSol );
	
}

void testSample2()
{
	MatrixDouble testMatrix = readMatD( "samplematrix2.txt", 4, 4 );
	ArrayDouble colVec = buildArrD(4, -4.0, 1.0, -3.0, 5.0);
	ArrayDouble testSol = forwSubst( testMatrix, colVec );
	
	printArrDPar( testSol );
	
	freeArrD(colVec);
	freeMatD( testMatrix );
	freeArrD( testSol );
	
}

void testSample1Back( )
{
	MatrixDouble testMatrix = readMatD( "samplematrix1back.txt", 3, 3 );
	ArrayDouble colVec = buildArrD(3, 1.0, 1.0, 6.0);
	ArrayDouble testSol = backSubst( testMatrix, colVec );
	
	printArrDPar( testSol );
	
	freeArrD(colVec);
	freeMatD( testMatrix );
	freeArrD( testSol );
	
}

void testSample2Back()
{
	MatrixDouble testMatrix = readMatD( "samplematrix2back.txt", 4, 4 );
	ArrayDouble colVec = buildArrD(4, 4.0, 1.0, 1.0, 5.0);
	ArrayDouble testSol = backSubst( testMatrix, colVec );
	
	printArrDPar( testSol );
	
	freeArrD(colVec);
	freeMatD( testMatrix );
	freeArrD( testSol );
	
}

int main()
{
	testSample2Back();
	
	
	return 0;
}
