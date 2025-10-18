// Es5_1.c LU Decomposition

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "DoubleUtilities.c"

#define PI 3.1415926535897932384

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

