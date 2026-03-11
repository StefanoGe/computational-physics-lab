#ifndef LINEARSYS
#define LINEARSYS

#include "matrix.h"
#include "array.h"
#include "genutil.h"

int linst_forwsubst_inplace(const Matrix *L, Array *b, bool is_diag_one);

int linst_forwsubst(const Matrix *L, const Array *b, Array *sol, bool is_diag_one);

int linst_backsubst_inplace( const Matrix *U, const Array *b);

int linst_backsubst( const Matrix *U, const Array *b, Array *sol);

int linst_lup_factor(Matrix *A, int *pivots, double tol);

double linst_lup_det(const Matrix *lu, int *pivots);

int linst_lup_solve_inplace(const Matrix *lu, int *pivots, Array *b);

Matrix linst_lu_extract_u(const Matrix *lu);

Matrix linst_lu_extract_l(const Matrix *lu);

Matrix linst_lu_extract_p(int *pivots, int dim);

int linst_lu_factor_no_pivot(Matrix *A, double tol);

void linst_lsqr_lup ( const Matrix *A, const Array *b, Array *x );

typedef struct {
    ParamFunc *funcs;   
    int nfuncs;          
} LinearModel;

LinearModel linear_model_new(ParamFunc *funcs, int nfuncs);

LinearModel linear_model_alloc(int nbasis);

void linear_model_free(LinearModel *m);

void linst_lsqr_fit_linear( const Array *x, const Array *y, 
const LinearModel *model, Array *coeffs );

#endif
