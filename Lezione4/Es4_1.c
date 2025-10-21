#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "DoubleUtilities.c"
#include "LinearSystems.c"

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
