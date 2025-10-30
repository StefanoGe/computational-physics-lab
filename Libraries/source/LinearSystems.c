// LinearSystems.c

#ifndef SG_LINEAR_SYSTEMS_C
#define SG_LINEAR_SYSTEMS_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comp_physics.h"

#define _FORW_SUBST 0
#define _BACK_SUBST 1


Linear_Model alloc_linear_model( int npar )
{
	Linear_Model model;
	model.npar = npar;
	model.funcs = ( Func_Ptr * ) malloc( sizeof(Func_Ptr) * npar );
	
	return model;
}

static inline void checkMatrixColHeight(MatrixDouble matrix, ArrayDouble colVec,
int func_caller)
{
	char * func = NULL;
	if( func_caller == _FORW_SUBST )
		func = "backSubst";
	else if( func_caller == _BACK_SUBST )
		func = "forwSubst";
	
	if( matrix.nrows != colVec.length )
		raiseErr( "Matrix and colVec should have same height!"
		"err in func %s.\nmatrix.nrows = %d, colVec.length = %d", func,
		matrix.nrows, colVec.length );
	
	if( matrix.nrows != colVec.length )
		raiseErr( "Matrix should be square!"
		"err in func %s.\nmatrix.nrows = %d, matrix.ncols = %d", func,
		matrix.nrows, matrix.ncols );
}

ArrayDouble forwSubst( MatrixDouble matrix, ArrayDouble colVec)
{
	checkMatrixColHeight( matrix, colVec, _FORW_SUBST );
	
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

MatrixDouble forw_sub_mat( MatrixDouble L, MatrixDouble known_terms)
{
	const int dim = L.nrows;
	const int n_sol = known_terms.ncols;
	
	MatrixDouble solutions = allocMatD( dim, n_sol );
	
	double partSum = 0;
	
	for( int p = 0; p < n_sol; p ++ )
	{
		for( int row = 0; row < dim; row++ )
		{
			partSum = 0;
			for( int col = 0; col < row; col++ )
				partSum += solutions.val[col][p] * L.val[row][col];
			solutions.val[row][p] = 
				(known_terms.val[row][p] - partSum)/L.val[row][row];
		}
	}
	
	return solutions;
}

ArrayDouble backSubst( MatrixDouble matrix, ArrayDouble colVec)
{
	checkMatrixColHeight( matrix, colVec, _BACK_SUBST );
	
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

MatrixDouble back_sub_mat( MatrixDouble U, MatrixDouble known_terms)
{
	const int dim = U.nrows;
	const int n_sol = known_terms.ncols;
	
	MatrixDouble solutions = allocMatD( dim, n_sol );
	
	double partSum = 0;
	
	for( int p = 0; p < n_sol; p ++ )
	{
		for( int row = dim - 1; row >= 0; row-- )
		{
			partSum = 0;
			for( int col = row + 1; col < dim; col++ )
				partSum += solutions.val[col][p] * U.val[row][col];
			solutions.val[row][p] = 
				(known_terms.val[row][p] - partSum)/U.val[row][row];
		}
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

LUP_Mats alloc_LUP( int dim )
{
	LUP_Mats lumats;
	lumats.dim = dim;
	lumats.L = allocMatD(dim, dim);
	lumats.U = allocMatD(dim, dim);
	lumats.P = allocMatD(dim, dim);
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

static inline void LUD_helper( LUMats lumats, MatrixDouble Acp, int n, int dim )
{
	// Build U
	for( int col = n; col < dim; col++ )
		lumats.U.val[n][col] = Acp.val[n][col];
	
	// Buil L
	for( int row = n + 1; row < dim; row ++ )
		lumats.L.val[row][n] = Acp.val[row][n] / lumats.U.val[n][n];
	
	// Subtract over Acp
	for( int row = 0; row < dim; row ++ )
		for( int col = 0; col < dim; col ++ )
	// Col/lines with index less/equal than n are set to 0 
			Acp.val[row][col] -= ( row <= n || col <= n ) ? 0 :
				lumats.L.val[row][n] * lumats.U.val[n][col];
}

LUMats LUDecompLow (MatrixDouble A)
{
	if( A.nrows != A.ncols )
		raiseErr("A must be square! In func LUDecomp\n"
		"A.nrows = %d, A.ncols = %d", A.nrows, A.ncols);
	
	const int dim = A.nrows;
	
	LUMats lumats = allocLUMats( dim );
	MatrixDouble Acp = allocMatD( dim, dim );
	cpMatToMatD( A, Acp, 0, 0 );
	
	//printMatDGraph( Acp );
	
	diagMatD(lumats.U, 0);
	diagMatD(lumats.L, 1);
	
	for (int n = 0; n < dim; n ++)
		LUD_helper( lumats, Acp, n, dim );
	
	freeMatD( Acp );
	
	return lumats;
}

inline static int find_pivot(MatrixDouble A_cp, int n)
{
	double max = fabs(A_cp.val[0][n]);
	double current;
	int max_index = 0;
	
	for( int row = 1; row < A_cp.nrows; row++ )
		if( max <= ( current = fabs(A_cp.val[row][n] ) ) )
		{	
			max = current;
			max_index = row;
		}
		
	printMatDGraph(A_cp);
	
	fprintf(stderr,"max_index = %d\n", max_index);
	
	if( fabs(max) <= 1e-6 )
		return -1;
	
	return max_index;
}

inline static void LUP_D_permuter ( LUP_Mats mats, MatrixDouble A_cp, int n )
{
	int pivot_row = find_pivot(A_cp, n);
	
	if( pivot_row == -1 )
		raiseErr( "Matrix of which solution is to be found is singular!\n" );
	
	if( pivot_row == n )
		return;
		
	exchange_rows( A_cp, n, pivot_row );
	//exchange_rows( mats.L, n, pivot_row );
	exchange_rows( mats.P, n, pivot_row );
	
}

inline static void LUP_D_helper_bad( LUP_Mats mats, MatrixDouble Acp, int n, int dim )
{
	LUP_D_permuter(mats, Acp, n);
	
	// Build U
	for( int col = n; col < dim; col++ )
		mats.U.val[n][col] = Acp.val[n][col];
	
	// Build L
	for( int row = n + 1; row < dim; row ++ )
		mats.L.val[row][n] = Acp.val[row][n] / mats.U.val[n][n];
	
	// Subtract over Acp
	for( int row = 0; row < dim; row ++ )
		for( int col = 0; col < dim; col ++ )
	// Col/lines with index less/equal than n are set to 0
			Acp.val[row][col] -= mats.L.val[row][n] * mats.U.val[n][col];
}

inline static void LUP_D_helper( LUP_Mats mats, MatrixDouble Acp, int n, int dim )
{
	int pivot = find_pivot( Acp, n );
	
	if( pivot == -1 )
		raiseErr( "Matrix of which solution is to be found is singular!\n" );
	
	exchange_rows( mats.P, n, pivot );
	
	// Build U
	for( int col = n; col < dim; col++ )
		mats.U.val[n][col] = Acp.val[pivot][col];
	
	// Build L
	for( int row = n + 1; row < dim; row ++ )
		mats.L.val[row][n] = Acp.val[row][n] / mats.U.val[n][n];
		
	// Subtract on A
		
	for( int row = 0; row < dim; row++ )
		for( int col = 0; col < dim; col++ )
			Acp.val[row][col] -= mats.L.val[row][n] * mats.U.val[n][col];
}


LUP_Mats LUP_decomposition ( MatrixDouble A )
{
	if( A.nrows != A.ncols )
		raiseErr("A must be square! In func LUDecomp\n"
		"A.nrows = %d, A.ncols = %d", A.nrows, A.ncols);
	
	const int dim = A.nrows;
	
	LUP_Mats lup_mats = alloc_LUP( dim );
	MatrixDouble A_cp = allocMatD( dim, dim );
	cpMatToMatD( A, A_cp, 0, 0 );
	
	//printMatDGraph( A_cp );
	
	diagMatD(lup_mats.U, 0);
	diagMatD(lup_mats.L, 1);
	diagMatD(lup_mats.P, 1);
	
	for (int n = 0; n < dim; n ++)
		LUP_D_helper( lup_mats, A_cp, n, dim );

	matMultD( lup_mats.P, lup_mats.L, &lup_mats.L );
	
	freeMatD( A_cp );
	
	return lup_mats;
}

ArrayDouble solve_LU( MatrixDouble A, ArrayDouble known_terms )
{
	LUMats lumat = LUDecompLow(A);
	
	ArrayDouble y = forwSubst( lumat.L, known_terms );
	ArrayDouble x = backSubst( lumat.U, y );
	
	freeAllMatD( lumat.L, lumat.U, NULL_MAT );
	freeArrD( y );
	
	return x;
}

double detTri(MatrixDouble triMat)
{
	double determinant = 1;
	for(int i = 0; i < triMat.nrows; i++)
		determinant *= triMat.val[i][i];
	
	return determinant;
}


double vec_1norm( ArrayDouble vec )
{
	double sum = 0;
	for( int i = 0; i < vec.length; i ++ )
		sum+= fabs(vec.val[i]);
	
	return sum;
}

double vec_infnorm( ArrayDouble vec )
{
	double max = fabs(vec.val[0]);
	double current;
	
	for(int i = 1; i < vec.length; i ++)
		if( max < ( current = vec.val[i] ) )
			max = current;
	
	return max;
}

double mat_infnorm( MatrixDouble A )
{
	double max = 0;
	double current;

	for( int i = 0; i < A.nrows; i ++ )
	{
		current = vec_1norm( getRowMatD( A, i ) );
		if (current > max)
			max = current;
	}
	
	return max;
}

double mat_1norm( MatrixDouble A )
{
	double max = 0;
	double sum = 0;

	for( int j = 0; j < A.ncols; j ++ )
	{
		sum = 0;
		for( int i = 0; i < A.nrows; i++  )
			sum += A.val[i][j];
		
		if( max < sum )
			max = sum;	
	}
	
	return max;
}

MatrixDouble mat_inv ( MatrixDouble A )
{
	
	LUMats lumats = LUDecompLow( A );
	MatrixDouble id = allocMatD( A.nrows, A.ncols );
	diagMatD( id, 1 );
	
	eprint("ok 1\n");
	
	MatrixDouble inv_U = back_sub_mat( lumats.U, id );
	MatrixDouble inv_L = forw_sub_mat( lumats.L, id );
	
	eprint("ok 2\n");
	
	MatrixDouble inv = matMultD( inv_U, inv_L, CREATE_MAT );
	
	eprint("ok 3\n");
	
	freeAllMatD( lumats.U, lumats.L, id, inv_L, inv_U, NULL_MAT );
	
	eprint("ok 3.5\n");
	
	return inv;
}

ArrayDouble min_x ( MatrixDouble A, ArrayDouble b )
{
	//We should have more equations than variables
	if( A.nrows < A.ncols )
		raiseErr( "A.nrows = %d should be less than A.ncols = %d\n", A.nrows, A.ncols );
	
	//Number of equations should match height of the known term vector
	if( A.nrows != b.length )
		raiseErr( "A.nrows = %d and b.length = %d should have same length.\n", A.nrows, b.length );
	
	MatrixDouble Atr = transposeMatD( A, CREATE_MAT );
	MatrixDouble Atr_A = matMultD( Atr, A, CREATE_MAT );
	ArrayDouble Atr_b = mat_vec_mult( Atr, b, NULL );
	LUMats lumat = LUDecompLow( Atr_A );
	ArrayDouble z = forwSubst( lumat.L, Atr_b );
	ArrayDouble x = backSubst( lumat.U, z );
	
	
	freeAllMatD( Atr, Atr_A, lumat.U, lumat.L, NULL_MAT );
	freeAllArrD( z, Atr_b, NULL_ARR );
	
	return x;
}

ArrayDouble linear_ls_pol_fitting( 
							ArrayDouble x_data, ArrayDouble y_data, int max_n )
{
// Function for the polynomial fitting of data by the least squares method. 
// max_n fixes the highest degree of the polynomials used.

	if( x_data.length != y_data.length )
		raiseErr( "x_data = %d and y_data = %d must have same size\n", 
			x_data.length, y_data.length );
			
	if( max_n < 0 )
		raiseErr(  );
			
	const int n_data = x_data.length;
	const int n_parameters = max_n +1;
	MatrixDouble A = allocMatD( n_data, n_parameters );
	
	// Set first col to 1 ( as it is the col of the constant function)
	for(int row = 0; row < n_data; row++ )
		A.val[row][0] = 1.0;
	
	for( int col = 1; col < n_parameters; col++ )
		for( int row = 0; row < n_data; row++ )
			A.val[ row ][ col ] = A.val[row][col - 1] * x_data.val[row];
			
	ArrayDouble result = min_x( A, y_data );
	
	freeMatD(A);
	
	return result;
}
//idea: creo un algoritmo che mi permette di scegliere se usare QR, LU or Q-less !!!
ArrayDouble linear_least_square_fitting( ArrayDouble x_data, ArrayDouble y_data,
			Linear_Model model)
{
	if( x_data.length != y_data.length )
		raiseErr( "must have same size\n" );
		
	const int n_parameters = model.npar;
	const int n_data = x_data.length;
		
	MatrixDouble A = allocMatD( n_data, n_parameters );
	
	//Costruisco la matrice di f
	
	for( int col = 0; col < n_parameters; col++ )
		for( int row = 0; row < n_data; row++ )
			A.val[ row ][ col ] = (*model.funcs[ col ])(x_data.val[row]);
	
	//La do in pasto a min_x
	
	ArrayDouble result = min_x( A, y_data );
	
	freeMatD(A);
	
	//Ho il risultato
	
	return result;
}

Param_Linear_Model alloc_param_linear_model( int nfuncs ){
	
	Param_Linear_Model model;
	model.funcs =  malloc( nfuncs * sizeof( Par_Func ) );
	model.nfuncs = nfuncs;
	
	return model;
}


QR_Mats alloc_QR( int nrows, int ncols  )
{
	QR_Mats qr_mats;
	qr_mats.Q = allocMatD( nrows, ncols );
	qr_mats.R = allocMatD( ncols, ncols );
	return qr_mats;
}

static inline double norm( MatrixDouble * Q, int col )
{
	double norm_sq = 0;
	for( int i = 0; i < Q->nrows; i++ )
		norm_sq += Q->val[i][col] * Q->val[i][col];
		
	return sqrt(norm_sq);
}

static inline double normalise( MatrixDouble * Q, int col )
{		
	double vec_norm = norm( Q, col );
	for( int i = 0; i < Q->nrows; i++ )
		Q->val[i][col] /= vec_norm;
		
	return vec_norm;
}

static inline void QR_decomp_helper( const QR_Mats * qr, const int i, const int j,
const int height )
{
	double scalar_product = 0;
	for( int p = 0; p < height; p++ )
		scalar_product += qr->Q.val[p][i] * qr->Q.val[p][j];
	
	qr->R.val[i][j] = scalar_product;
	
	for( int p = 0; p < height; p++ )
	{
		qr->Q.val[p][j] -= scalar_product * qr->Q.val[p][i];
	}
}

// Sistemare meglio questa funzione

QR_Mats QR_decomp( MatrixDouble A )
{
	const int height = A.nrows;
	const int nvecs = A.ncols;
	QR_Mats qr = alloc_QR( height, nvecs );
	diagMatD(qr.R, 0);
	cpMatToMatD( A, qr.Q, 0, 0 );
	
	// loop che controlla il q che sto usando
	for( int i = 0; i < nvecs; i ++ )
	{
		qr.R.val[i][ i ] = normalise( &qr.Q, i );
		// Mi serve ora un algoritmo che calcoli i vari vj
		for( int j = i + 1; j < nvecs; j++)
		{
			QR_decomp_helper( &qr, i, j, height );
		}
	}
	return qr;
}

#endif











