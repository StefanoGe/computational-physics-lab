// Es4Opt.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "DoubleUtilities.c"

void checkMatrixsHeight(MatrixDouble matrix, MatrixDouble knownTerms)
{
	if( matrix.nrows != knownTerms.nrows )
	{
		printf("Matrix and knownTerms should have same height! err in func forwSubst.\n");
		exit(EXIT_FAILURE);
	}
	if( matrix.ncols != matrix.nrows )
	{
		printf("La matrice deve essere quadrata! err in func forwSubst.\n");
		exit(EXIT_FAILURE);
	}
}

ArrayDouble forwSubstGen( MatrixDouble matrix, ArrayDouble knownTerms)
{
	checkMatrixsHeight( matrix, knownTerms );
	
	const int n = knownTerms.nrows;
	const int p = knownTerms.ncols;
	
	MatrixDouble solutions = allocMatD( n, p );
	
	ArrayDouble partSums = allocArrD( p );
	
	for( int row = 0; row < n; row++ )
	{
		setValueArrD( partSums, 0 );
		for( int col = 0; col < row; col++ )
		{
			carrTimesScalarD(solutions.val[col],p, matrix.val[row][col]);
			ArrayDPlusCarr
		}
		solutions.val[row] = (colVec.val[row] - partSum)/matrix.val[row][row];
	}
	
	return solutions;
}


int main()
{
	
	
	return 0;
}
