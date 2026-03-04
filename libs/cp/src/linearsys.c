#include "linearsys.h"
#include "genutil.h"
#include <math.h>

/*
 * Forward substitution linear system solver.
 * Finds unknown column vector of equation LX=b, with L lower triangular matrix.
 * 
 * args:
 * L - lower triangular matrix
 * b - right hand side
 * sol - solution destination
 */

int linst_forwsubst(const Matrix *L, const Array *b, Array *sol, bool is_diag_one)
{
	if(!mat_is_squared(L))
		raiseErr("L must be squared");
	
	if( sol->size != L->ncols )
		arr_init(sol, L->ncols);
	
	arr_cp(b, sol);
	
	return linst_forwsubst_inplace(L, sol, is_diag_one);
}

int linst_forwsubst_inplace(const Matrix *L, Array *b, bool is_diag_one)
{	
	if(L->ncols!=L->nrows)
		raiseErr("L must be squared");
		
	const int dim = L->ncols;
	
	if(b->size != dim)
		raiseErr("L and b must have same dimension");
	
	double partial_sum;
	
	for( int row=0; row < dim; row++ )
	{
		if (MATP(L,row,row) == 0.0)
			return row;
		
		partial_sum = 0;
		for( int col = 0; col < row; col++ )
			partial_sum += ARRP(b,col) * MATP(L,row,col);
		ARRP(b,row) = (ARRP(b,row) - partial_sum);
		if(!is_diag_one)
			ARRP(b,row)/=MATP(L,row,row);
	}
	return 0;
}

/*
 * Backward substitution linear system solver.
 * Finds unknown column vector of equation UX=b, with U upper triangular matrix.
 * 
 * args:
 * U - upper triangular matrix
 * b - right hand side
 * sol - solution destination
 */

int linst_backsubst( const Matrix *U, const Array *b, Array *sol)
{
	if(!mat_is_squared(U))
		raiseErr("U must be squared");
	
	if( sol->size != U->ncols )
		arr_init(sol, U->ncols);
	
	arr_cp(b, sol);
	
	return linst_backsubst_inplace(U, sol);
}


int linst_backsubst_inplace( const Matrix *U, const Array *b)
{
	if(U->ncols!=U->nrows)
		raiseErr("U must be squared");
		
	const int dim = U->ncols;
	
	if(b->size != dim)
		raiseErr("U and b must have same dimension");
	
	double partial_sum;
	
	for( int row = dim - 1; row >= 0; row-- )
	{
		if (MATP(U,row,row) == 0.0)
			return row;
		
		partial_sum = 0;
		for( int col = row + 1; col < dim; col++ )
			partial_sum += ARRP(b,col) * MATP(U,row,col);
		ARRP(b,row) = (ARRP(b,row) - partial_sum)/MATP(U,row,row);
	}
	return 0;
}

LUP_Result lup_init(int size)
{
	LUP_Result result={0};
	result.size=size;
	result.L=mat_new(size, size);
	result.U=mat_new(size,size);
	mat_diag(&result.U, 0);
	mat_diag(&result.L, 1);
	return result;
}

static inline int find_pivot(const Matrix *A, int n, int dim)
{
	double max_pivot=0;
	int pivot_row=n;
	for(int row=n; row<dim; row++)
	{
		if(fabs(MATP(A,row,n))>=max_pivot)
		{
			max_pivot=fabs(MATP(A,row,n));
			pivot_row=row;
		}
	}
	return pivot_row;
}

/*
 * In place LU-factorization with row-pivoting
 * Factorization happens in-place, i.e. output is same as input
 * A. Diagonal ones of L are omitted.
 * 
 * args:
 * A - matrix to be factored
 * pivots - integer array into which pivots position are stored
 * tol - tolerance for checking when pivots are 0. if tol is set <0, 
 * 		default tolerance of 1e-14 is used instead.
 * 
 * returns:
 * 0 on success, n when n-th |pivot|<tol
 */

int linst_lup_factor(Matrix *A, int *pivots, double tol)
{	
	if(!mat_is_squared(A))
		raiseErr("A must be squared");
		
	if(tol<0)
		tol=1e-14;
	
	const int dim = A->ncols;
	double pivot = 0;
	
	for(int n = 0; n < dim; n++)
	{
		// Find pivot
		pivots[n]=find_pivot(A, n, dim);
		if(pivots[n]!=n)
			mat_swap_row(A,n,pivots[n]);
		pivot=MATP(A,n,n);
		// Check if pivot is 0
		if(fabs(pivot)<=tol)
			return n;
		// Write L
		for(int row=n+1; row<dim; row++)
			MATP(A,row,n)=MATP(A,row,n)/pivot;
		// Subtract in A
		for(int row=n+1; row<dim; row++)
			for(int col=n+1; col<dim; col++)
				MATP(A, row, col) = MATP(A,row,col)-MATP(A,row,n)*MATP(A,n,col);
	}
	return 0;
}

double linst_lup_det(const Matrix *lu, int *pivots)
{
	double det=1;
	for(int i=0; i<lu->nrows; i++)
		det*=MATP(lu,i,i)*(pivots[i]==i ? 1 : -1);
	return det;
}

int linst_lup_solve_inplace(const Matrix *lu, int *pivots, Array *b)
{
	for(int i=0; i<b->size; i++)
		if(i!=pivots[i])
			SWAP( ARRP(b, i), ARRP(b,pivots[i]) );
	int info=0;
	
	if( (info=linst_forwsubst_inplace(lu, b, true) )>0)
		return -info;
	return linst_backsubst_inplace(lu,b);
}

Matrix linst_lu_extract_u(const Matrix *lu)
{
	Matrix u = mat_new(lu->nrows, lu->ncols);
	for(int i=0; i<lu->nrows; i++)
		for(int j=0; j<lu->ncols; j++)
			MAT(u,i,j)=( j>=i ? MATP(lu,i,j) : 0 );
	return u;
}

Matrix linst_lu_extract_l(const Matrix *lu)
{
	Matrix l = mat_new(lu->nrows, lu->ncols);
	for(int i=0; i<lu->nrows; i++)
		for(int j=0; j<lu->ncols; j++)
			MAT(l,i,j)=( j<i ? MATP(lu,i,j) : (j==i ? 1 : 0) );
	return l;
}

Matrix linst_lu_extract_p(int *pivots, int dim)
{
	Matrix P = mat_new(dim, dim);
	mat_diag(&P, 1);
	for(int i=0; i<dim; i++)
		if(pivots[i]!=i) mat_swap_row(&P, i, pivots[i]);
	
	return P;
}


