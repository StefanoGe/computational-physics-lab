// LinearSystems.c

#ifndef SG_LINEAR_SYSTEMS_C
#define SG_LINEAR_SYSTEMS_C

typedef struct {
	
	MatrixDouble U;
	MatrixDouble L;
	int dim;
	
} LUMats;

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

LUMats allocLUMats( int dim )
{
	LUMats lumats;
	lumats.dim = dim;
	lumats.L = allocMatD(dim, dim);
	lumats.U = allocMatD(dim, dim);
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
	
	//dividMat = allocMatD( dim, dim );
	
	printf("Upper is: \n");
	printMatDGraph( lumats.U );
	printf("Lower is: \n");
	printMatDGraph( lumats.L );
	
	return lumats;
}

LUMats LUDecompLow (MatrixDouble A)
{
	if( A.nrows != A.ncols )
		raiseErr("A must be square! In func LUDecomp\n");
	
	const int dim = A.nrows;
	
	LUMats lumats = allocLUMats( dim );
	MatrixDouble Acp = allocMatD( dim, dim );
	cpMatToMatD( A, Acp, 0, 0 );
	
	//printMatDGraph( Acp );
	
	diagMatD(lumats.U, 0);
	diagMatD(lumats.L, 1);
	
	for (int n = 0; n < dim; n ++)
	{
		// Costruisco U
		for( int i = n; i < dim; i++ )
			lumats.U.val[n][i] = Acp.val[n][i];
		
		// Costruisco L
		for( int j = n + 1; j < dim; j ++ )
			lumats.L.val[j][n] = Acp.val[j][n] / lumats.U.val[n][n];
			
		// Sottraggo su Acp
		for( int i = 0; i < dim; i ++ )
			for( int j = 0; j < dim; j ++ )
				Acp.val[i][j] -= lumats.L.val[i][n] * lumats.U.val[n][j];
	}
	
	freeMatD( Acp );
	
	return lumats;
}

#endif
