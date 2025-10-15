// Es5_1.c LU Decomposition

#include <stdio.h>
#include <stdlib.h>
#include "DoubleUtilities.c"

typedef struct {
	
	MatrixDouble U;
	MatrixDouble L;
	int dim;
	
} LUMats;

LUMats allocLUMats( int dim )
{
	LUMats lumats;
	lumats.dim = dim;
	lumats.L = allocMatD(dim, dim);
	lumats.U = allocMatD(dim, dim);
	setValueMatD(lumats.L, 0);
	setValueMatD(lumats.U, 0);
	return lumats;
}


LUMats LUDecomp( MatrixDouble A )
{
	if( A.nrows != A.ncols )
		raiseErr("A must be square! In func LUDecomp\n");
		
	const int dim = A.nrows;
		
	LUMats lumats = allocLUMats( dim );
		
	// First row of U is just first row of A
	
	cpArrayToRowMatD( lumats.U, 0, getRowMatD(A, 0) );
	
	// First col of L is first col of A / u11
	
	dividMat = allocMatD( dim, dim );
	divMat(  )
	
	printf("Upper is: \n");
	printMatDGraph( lumats.U );
	printf("Lower is: \n");
	printMatDGraph( lumats.L );
	
	return lumats;
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

int main()
{
	printf("Hello world!");
	
	MatrixDouble testMat = readMatD("samplematrix1.txt", 3, 3);
	LUMats lumats = LUDecomp ( testMat );
	
	
	freeMatD(testMat);
	freeMatD(lumats.U);
	freeMatD(lumats.L);
	
	return 0;
}

