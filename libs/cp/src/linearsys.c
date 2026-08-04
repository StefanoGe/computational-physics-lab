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


int linst_backsubst_inplace(const Matrix *U, Array *b)
{
	if(U->ncols!=U->nrows)
		raiseErr("U must be squared");
		
	const int dim = U->ncols;
	
	if(b->size != dim)
		raiseErr("U and b must have same dimension");
	
	return linst_backsubst_inplace_n(U,b,dim);
}


int linst_backsubst_inplace_n( const Matrix *U, Array *b, int dim)
{
	if(U->nrows<dim||U->ncols < dim)
		raiseErr("Matrix should be bigger nrows %d ncols %d dim %d",
			U->nrows, U->ncols, dim);
	
	if(b->size < dim)
		raiseErr("b should be bigger size %d dim %d", b->size, dim);
	
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

int linst_backsubst_bulk(const Matrix *U, const Matrix *B, Matrix *X)
{
	if(U->ncols!=U->nrows)
		raiseErr("U must be squared");
	
	const int dim = U->nrows;
	
	if(B->nrows != dim)
		raiseErr("B nrows must match linsys dimension");

	const int nrhs = B->ncols;
	
	mat_init(X,dim,nrhs);
	
	Array temp_arr = arr_new(dim);
	
	for(int i=0; i<nrhs; i++)
	{
		mat_col_to_arr(B,&temp_arr,i);
		
		if(linst_backsubst_inplace(U,&temp_arr))
			return 1;
		
		for(int j=0; j<dim; j++)
			MATP(X,j,i)=ARR(temp_arr,j);
	}
	
	arr_free(&temp_arr);
	
	return 0;
}

static inline int find_pivot(const Matrix *A, int n, int dim)
{
	double max_pivot=0;
	int pivot_row=n;
	for(int row=n; row<dim; row++)
		if(fabs(MATP(A,row,n))>max_pivot)
		{
			max_pivot=fabs(MATP(A,row,n));
			pivot_row=row;
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

int linst_lu_factor_no_pivot(Matrix *A, double tol)
{	
	if(!mat_is_squared(A))
		raiseErr("A must be squared");
		
	if(tol<0)
		tol=1e-14;
	
	const int dim = A->ncols;
	
	for(int n = 0; n < dim; n++)
	{
		// Check if pivot is 0
		if(fabs(MATP(A,n,n))<=tol)
			return n;
		// Write L
		for(int row=n+1; row<dim; row++)
			MATP(A,row,n)=MATP(A,row,n)/MATP(A,n,n);
		// Subtract in A
		for(int row=n+1; row<dim; row++)
			for(int col=n+1; col<dim; col++)
				MATP(A, row, col) = MATP(A,row,col)-MATP(A,row,n)*MATP(A,n,col);
	}
	return 0;
}

int linst_lup_solve_inplace(const Matrix *lu, const int *pivots, Array *b)
{
	for(int i=0; i<b->size; i++)
	{
//		eprint("pivots[i] = %d",pivots[i]);
		if(i!=pivots[i])
			SWAP( ARRP(b, i), ARRP(b,pivots[i]) );

	}
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

static inline void mat_atb(const Matrix *A, const Array *b, Array *x)
{
	const int dim_g = A->ncols;
	for (int i=0; i<dim_g; i++)
	{
		double sum = 0;

		for (int p=0; p<A->nrows; p++)
			sum += MATP(A,p,i) * ARRP(b,p);

		ARRP(x,i) = sum;
	}
}

void linst_lsqr_lup ( const Matrix *A, const Array *b, Array *x )
{
	//We should have more equations than variables
	if( A->nrows < A->ncols )
		raiseErr( "A.nrows = %d should be greater than or equal A.ncols = %d\n", 
			A->nrows, A->ncols );
	
	//Number of equations should match height of the right hand side vector
	if( A->nrows != b->size )
		raiseErr( "A.nrows = %d and b.length = %d should have same length.\n", 
			A->nrows, b->size );
	
	const int dim_g = A->ncols;
	Matrix G = mat_new(dim_g, dim_g); // Gram matrix A_tr A
	arr_init(x,dim_g);
	mat_gram(A,&G);
	mat_atb(A,b,x);
	
	
	int *pivots = malloc(sizeof(int)*dim_g);
	
	if( linst_lup_factor(&G, pivots, 0) )
		raiseErr("could not lup factor");
	
	if( linst_lup_solve_inplace(&G, pivots, x) )
		raiseErr("Could not solve linear system");
	
	
	free(pivots);
	mat_free(&G);
}

LinearModel linear_model_new(ParamFunc *funcs, int nfuncs) {
    LinearModel model;
    model.funcs  = funcs;
    model.nfuncs = nfuncs;
    return model;
}

LinearModel linear_model_alloc(int nfuncs) {
    LinearModel m;
    m.funcs = malloc(sizeof(ParamFunc) * nfuncs);
    m.nfuncs = nfuncs;
    return m;
}

void linear_model_free(LinearModel *m) {
    free(m->funcs);
    m->funcs=nullptr;
    m->nfuncs=0;
}

static inline void make_linear_model_matrix(Matrix *A, const LinearModel *model,
	const Array *x)
{
	
	const int n_parameters = model->nfuncs;
	const int n_data = x->size;
	
	mat_init(A,n_data, n_parameters );
	
	for( int col = 0; col < n_parameters; col++ )
		for( int row = 0; row < n_data; row++ )
			MATP(A,row,col) = eval(model->funcs+col,ARRP(x,row));

	
}

void linst_lsqr_fit_linear_lup( const Array *x, const Array *y, 
	const LinearModel *model, Array *coeffs )
{
	if( x->size != y->size )
		raiseErr( "must have same size\n" );
	
	Matrix A={0};
	
	make_linear_model_matrix(&A, model, x);
	
	// Apply least square solver
	
	linst_lsqr_lup( &A, y, coeffs );
	
	mat_free(&A);
}

double linst_linear_eval(double x, const LinearModel *l, const Array *coeffs)
{
	double sum = 0;
	for(int i=0; i<l->nfuncs; i++)
		sum+= ARRP(coeffs, i) * eval(l->funcs+i, x);
	return sum;
}



void linst_qr_mgs(const Matrix *A, Matrix *Q, Matrix *R)
{
	const int vsize=A->nrows;
	const int nvecs=A->ncols;
	mat_init(Q,vsize,nvecs);
	mat_init(R,nvecs,nvecs);
	
	mat_cp(A,Q);
	mat_diag(R,0);
	
	for(int s=0; s<nvecs; s++)
	{
		double qi_norm = 0;
		for(int p=0; p<vsize; p++)
			qi_norm += MATP(Q,p,s)*MATP(Q,p,s);
		qi_norm=sqrt(qi_norm);
        MATP(R,s,s) = qi_norm;
		for(int p=0; p<vsize; p++)
			MATP(Q,p,s)/=qi_norm;
		
		for(int v=s+1; v<nvecs; v++)
		{
			double dot = 0;
			for(int p = 0; p<vsize; p++)
				dot+=MATP(Q,p,s)*MATP(Q,p,v);	
			MATP(R,s,v)=dot;
			
			for(int p=0; p<vsize; p++)
				MATP(Q,p,v)-=dot*MATP(Q,p,s);
		}
		
		
	}
}

void linst_lsqr_qr( const Matrix *A, const Array *b, Array *x )
{
	//We should have more equations than variables
	if( A->nrows < A->ncols )
		raiseErr( "A.nrows = %d should be greater than or equal A.ncols = %d\n", 
			A->nrows, A->ncols );
	
	//Number of equations should match height of the right hand side vector
	if( A->nrows != b->size )
		raiseErr( "A.nrows = %d and b.length = %d should have same length.\n", 
			A->nrows, b->size );
	
	const int dim_g = A->ncols;
	Matrix Q={0};
	Matrix R={0};
	linst_qr_mgs(A,&Q,&R);

	arr_init(x,dim_g);
	mat_atb(&Q,b,x);
	
	linst_backsubst_inplace(&R, x);
	
	mat_free(&R);
	mat_free(&Q);
}

void linst_lsqr_fit_linear_qr( const Array *x, const Array *y, 
	const LinearModel *model, Array *coeffs )
{
	if( x->size != y->size )
		raiseErr( "must have same size\n" );
	
	Matrix A={0};
	
	make_linear_model_matrix(&A, model, x);
	
	// Apply least square solver
	
	linst_lsqr_qr( &A, y, coeffs );
	
	mat_free(&A);
}

void linst_qless(const Matrix *A, Matrix *Q, Matrix *R, const Array *b,
	Array *z)
{
	const int vsize=A->nrows;
	const int nvecs=A->ncols;
	mat_init(Q,vsize,nvecs);
	mat_init(R,nvecs,nvecs);
	arr_init(z,nvecs);
	
	mat_cp(A,Q);
	mat_diag(R,0);
	
    // working copy of b (will become z progressively)
    Array w={0};
    arr_init(&w, vsize);
    arr_cp(b, &w);
	
	for(int s=0; s<nvecs; s++)
	{
		double qi_norm = sqrt(mat_dot_col_col(Q,s,s));
		if(qi_norm <= 1e-15)
			raiseErr("Linearly dependent columns: norm=%g", qi_norm);
			
        MATP(R,s,s) = qi_norm;
		mat_col_scale(Q,s,1.0/qi_norm);
		
		double zs=mat_dot_col_arr(Q,s,&w);
		ARRP(z,s)=zs;
		
		for(int v=s+1; v<nvecs; v++)
		{
			double dot = mat_dot_col_col(Q,s,v);
			MATP(R,s,v)=dot;
			
			mat_axpy_col_col(Q,v,-dot,Q,s);
		}

        // augmented behavior
        for (int p = 0; p < vsize; p++)
            ARR(w,p) -= zs * MATP(Q,p,s);
	}
	arr_free(&w);
}
/*
void linst_lsqr_qless( const Matrix *A, const Array *b, Array *x )
{
	//We should have more equations than variables
	if( A->nrows < A->ncols )
		raiseErr( "A.nrows = %d should be greater than or equal A.ncols = %d\n", 
			A->nrows, A->ncols );
	
	//Number of equations should match height of the right hand side vector
	if( A->nrows != b->size )
		raiseErr( "A.nrows = %d and b.length = %d should have same length.\n", 
			A->nrows, b->size );
	
	
	Matrix Q={0};
	Matrix R={0};
	linst_qless(A,&Q,&R,b,x);
	
	linst_backsubst_inplace(&R, x);
	
	mat_free(&R);
	mat_free(&Q);
}
*/

void linst_lsqr_qless( const Matrix *A, const Array *b, Array *x )
{
    // build augmented matrix [A | b]
    Matrix Ab=mat_aug_arr_new(A,b);

    // QR decomposition
    Matrix R={0};
    Matrix Q={0};
    linst_qr_mgs(&Ab,&Q,&R);

	arr_init(x,A->ncols);

    for (int i = 0; i < A->ncols; i++)
        ARRP(x,i) = MATP(&R,i,A->ncols);

    // --- solve Rx = z ---
    if (linst_backsubst_inplace_n(&R, x,A->ncols))
        raiseErr("Back substitution failed");

    mat_free(&R);
    mat_free(&Ab);
    mat_free(&Q);
}

void linst_lsqr_fit_linear_qless( const Array *x, const Array *y, 
	const LinearModel *model, Array *coeffs )
{
	if( x->size != y->size )
		raiseErr( "must have same size\n" );
	
	Matrix A={0};
	
	make_linear_model_matrix(&A, model, x);
	
	// Apply least square solver
	
	linst_lsqr_qless( &A, y, coeffs );
	
	mat_free(&A);
}


/*
 * Pure QR algorithm to find the eigenvalues of a real symmetric matrix
 * whose eigenvalues have all different absolute values, and whose corresponding
 * eigenvector matrix has all nonsingular leading principal minors
 * 
 * args:
 * A: matrix satisfying the requirements above
 * eigenvalues: array into which to store the resulting eigenvalues
 * 
 */

double off_diagonal_norm(const Matrix *A)
{
	double off = 0;
	for (int i = 1; i < A->nrows; i++)
		for (int j = 0; j < i; j++)
			off += MATP(A,i,j) * MATP(A,i,j);

	off = sqrt(off);
	return off;
}

int linst_pure_qr(const Matrix *A, Array *eigenvalues, Matrix *evecs, 
	double tol)
{
	const int dim=A->nrows;
	
	Matrix Acp={0};
	mat_cp(A,&Acp);
	
	Matrix Q={0};
	Matrix R={0};
	mat_init(evecs,dim,dim);
	mat_diag(evecs,1);
	
	double off_norm=INFINITY;
	
	int i=0;
	
	while(tol<off_norm)
	{
		linst_qr_mgs(&Acp,&Q,&R);
		mat_mult(&R,&Q,&Acp);
		mat_mult(evecs,&Q,evecs);
		off_norm=off_diagonal_norm(&Acp);
		
		i++;
		if(i>=1000)
			raiseErr("loop did not stop");
	}
	
	arr_init(eigenvalues,dim);
	for(int p=0; p<dim; p++)
		ARRP(eigenvalues,p)=MAT(Acp,p,p);
	
	Matrix *to_free[]={&Q,&R,&Acp};
	mat_free_many(to_free,3);
	return i;
}

double linst_inf_norm(const Matrix *m)
{
	if(!mat_is_squared(m))
		raiseErr("mat must be squared");
	
	double max_sum = 0;
	double partial_sum=0;
	for(int row=0; row<m->nrows; row++)
	{
		partial_sum=0;
		for(int col=0; col<m->ncols;col++)
			partial_sum+=fabs(MATP(m,row,col));
		max_sum=MAX(partial_sum,max_sum);
	}
	
	return max_sum;
}

int linst_inv_qr(const Matrix *A, Matrix *invA)
{
	if(!mat_is_squared(A))
		raiseErr("mat must be squared");

	const int dim = A->nrows;

	mat_init(invA,dim,dim);

	Matrix Q={0};
	Matrix R={0};
	Matrix id=mat_new(dim,dim);
	
	mat_diag(&id,1);

	linst_qr_mgs(A,&Q,&R);
	
	mat_transpose(&Q,&Q);
	
	if(linst_backsubst_bulk(&R,&id,invA))
		return 1;
	
	mat_mult(invA, &Q, invA);
	
	Matrix *to_free[]={&Q,&R,&id};
	
	mat_free_many(to_free,3);
	
	return 0;
}

double linst_cond_numb_inf(const Matrix *A)
{
	if(!mat_is_squared(A))
		raiseErr("mat must be squared");
	
	double norm = linst_inf_norm(A);
	Matrix invA = mat_new(A->nrows,A->ncols);
	linst_inv_qr(A, &invA);
	
	double norm_inv=linst_inf_norm(&invA);
	
	return norm*norm_inv;
}


