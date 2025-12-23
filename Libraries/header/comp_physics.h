// My general header

#ifndef SG_COMP_PHYSICS_HEADER_H
#define SG_COMP_PHYSICS_HEADER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>

#define PI 3.1415926535897932384626433

#define CREATE_MAT NULL
#define NULL_ARR (ArrayDouble){NULL,0}
#define NULL_MAT (MatrixDouble){NULL, 0,0}
#define NULL_INFO (PlotInfo){NULL, NULL,0}

#ifndef GENUTIL_H
#define GENUTIL_H

#define UNUSED(expr) do { (void)(expr); } while(0)

#define raiseErr( msg, ...) \
	do { \
		fprintf(stderr, "[%s: line %d in %s]: " msg, __FILE__, __LINE__, __func__ \
		__VA_OPT__(,) __VA_ARGS__); \
		exit(EXIT_FAILURE); \
	} while (0)

#define eprint(a, ...) fprintf( stderr, a "\n" __VA_OPT__(,) __VA_ARGS__ )

FILE * openFile( const char * fileName, const char * mode );

#endif

typedef double(*Func_Ptr)(double);

typedef double(*ParamFuncPtr) (double, void *);

typedef struct {
	ParamFuncPtr param_func_ptr;
	void * params;
	int nparams;
} Par_Func;

double evaluate( const Par_Func * self, double x );

typedef struct _arrayInt
{
	int * val;
	int length;
} ArrayInt;

typedef struct {
	int x;
	int y;
} Point2DI;

ArrayInt allocArrI ( int length);

void freeArrI( ArrayInt array );

void printArrIPar (ArrayInt array);

typedef double(*Func_Ptr)(double);

typedef struct
{
	double * val;
	int length;
} ArrayDouble;

typedef struct MatrixDouble
{
	double ** val;
	int nrows;
	int ncols;
}MatrixDouble;

double pown( double x, int n );

ArrayDouble allocArrD ( int length);

void freeArrD( ArrayDouble array );

ArrayDouble buildArrD ( int length, ... );

ArrayDouble linspaceD ( double start, double end, int num_values );

ArrayDouble diffArrD( ArrayDouble array1, ArrayDouble array2);

void printArrDPar (ArrayDouble array, char * format);

void print_arr( ArrayDouble array, char * format, FILE * output );

void print_arr_tofile( ArrayDouble array, char * file_name);

void setValueArrD(ArrayDouble array, double value);

void cpArrD( ArrayDouble source, ArrayDouble destination );

MatrixDouble allocMatD ( int nrows, int ncols );

void printMatDGraph( MatrixDouble matrix );

MatrixDouble readMatD(char * filename, int nrows, int ncols);

void freeMatD( MatrixDouble matrix );

ArrayDouble getRowMatD( MatrixDouble matrix, int row );

void passArrayToRowMatD( MatrixDouble matrix, int row, ArrayDouble array );

void cpArrayToRowMatD( MatrixDouble matrix, int row, ArrayDouble array );

void cpArrayToColMatD( ArrayDouble array, MatrixDouble dest, int col,
		bool freeArr );

MatrixDouble transposeMatD( MatrixDouble matrix, MatrixDouble * dest);

void printDatMatD( MatrixDouble matrix, char * filename, char * format, 
	bool doTranspose);

void setValueMatD( MatrixDouble matrix, double value );

char * firstWS( char * s );

ArrayInt parseIndices( const char * select );

MatrixDouble sliceMatD( MatrixDouble matrix, char * rowSelect, char * colSelect );

void cpMatToMatD( MatrixDouble source, MatrixDouble dest, 
	int startRow, int startCol );

void divMat( MatrixDouble matrix, double num, MatrixDouble dest );

void diagMatD( MatrixDouble matrix, double value );

MatrixDouble matMultD( MatrixDouble mat1, MatrixDouble mat2, MatrixDouble * dest);

MatrixDouble asColMatrix( ArrayDouble arr, bool destroySource );

void freeAllArrD( ArrayDouble arr, ... );

# define free_arrD( ... ) freeAllArrD( __VA_ARGS__, NULL_ARR )

void freeAllMatD( MatrixDouble matrix, ... );

# define free_MatD( ... ) freeAllMatD( __VA_ARGS__, NULL_MAT )

void exchange_rows( MatrixDouble matrix, int row1, int row2 );

ArrayDouble mat_vec_mult( MatrixDouble matrix, ArrayDouble array, 
ArrayDouble * dest );

MatrixDouble mat_diffD( MatrixDouble mat1, MatrixDouble mat2, MatrixDouble * dest);

ArrayDouble build_from( ArrayDouble arr, Func_Ptr func );

typedef struct {
	
	MatrixDouble U;
	MatrixDouble L;
	int dim;
	
} LUMats;

typedef struct {
	
	MatrixDouble U;
	MatrixDouble L;
	MatrixDouble P;
	int dim;
	
} LUP_Mats;

typedef struct {

	Func_Ptr * funcs;
	int npar;

} Linear_Model;

typedef struct {
	
	Par_Func * funcs;
	int nfuncs;
	
} Param_Linear_Model;

Par_Func alloc_par_func( ParamFuncPtr, int );

Linear_Model alloc_linear_model( int npar );

Param_Linear_Model alloc_param_linear_model( int nfuncs );

ArrayDouble forwSubst( MatrixDouble matrix, ArrayDouble colVec);

MatrixDouble forw_sub_mat( MatrixDouble L, MatrixDouble known_terms);

ArrayDouble backSubst( MatrixDouble matrix, ArrayDouble colVec);

MatrixDouble back_sub_mat( MatrixDouble U, MatrixDouble known_terms);

LUMats allocLUMats( int dim );

LUP_Mats alloc_LUP( int dim );

LUMats LUDecomp( MatrixDouble A );

LUMats LUDecompLow (MatrixDouble A);

LUP_Mats LUP_decomposition ( MatrixDouble A );

double detTri(MatrixDouble triMat);

double vec_1norm( ArrayDouble vec );

double vec_infnorm( ArrayDouble vec );

double mat_infnorm( MatrixDouble A );

double mat_1norm( MatrixDouble A );

MatrixDouble mat_inv ( MatrixDouble A );

ArrayDouble min_x ( MatrixDouble A, ArrayDouble b );

ArrayDouble linear_ls_pol_fitting( 
							ArrayDouble x_data, ArrayDouble y_data, int max_n );

ArrayDouble linear_least_square_fitting( ArrayDouble x_data, ArrayDouble y_data,
			Linear_Model model);

typedef struct _arrayFloat
{
	float * val;
	int length;
	
} ArrayFloat;

ArrayFloat allocArrF ( int length);

ArrayFloat buildArrF ( int length, ... );

void freeArrF( ArrayFloat array );

void printArrFPar (ArrayFloat array);
void setValueArrF(ArrayFloat array, float value);


ArrayFloat diffArrF( ArrayFloat array1, ArrayFloat array2);

ArrayDouble solve_LU( MatrixDouble A, ArrayDouble known_terms );


typedef struct{
	MatrixDouble Q;
	MatrixDouble R;
} QR_Mats;

QR_Mats alloc_QR( int nrows, int ncols  );

QR_Mats QR_decomp( MatrixDouble A );

// vector.c

#define DEFAULT_VEC_SIZE 64

typedef struct {
	double * val;
	int size;
	int length;
}VectorD;

extern const VectorD NULL_VEC;

VectorD init_vecD();

VectorD alloc_vecD( int size );

void appendD( VectorD *, double );

void _free_all_vecD( VectorD *, ... );
# define free_vecD(vec, ...) \
	_free_all_vecD( vec, __VA_ARGS__ __VA_OPT__(,) &NULL_VEC );

VectorD build_vecD ( int length, ... );

void print_vecD( const VectorD *, char * format, FILE * output );

void std_print_vecD( const VectorD * );

VectorD vec_cp( const VectorD * source );

VectorD init_vec_length( int length );

VectorD vec_range( double x1, double x2, int num );

int vecD_necessary_size(int length);

VectorD vec_par_func( const VectorD *, const Par_Func * );

VectorD vec_full ( int length, double value );

double vec_min( const VectorD * );

double vec_max( const VectorD * );

// interp.c

typedef struct {
	VectorD points;
	VectorD weights;
} BaricFitter;

BaricFitter init_bar_fitter( const VectorD * points );

typedef struct {
	VectorD points;
	VectorD weights;
	VectorD f_values;
	Func_Ptr func;
} BarFit;

BaricFitter init_bar_fitter( const VectorD * points );

BarFit bar_fir( const BaricFitter *, Func_Ptr );

double barf_get_value( const BarFit * barf, double x );

BaricFitter bar_fitter_eq_init( double x1, double x2, int n_points );

VectorD cheb2_nodes_def( int num );

VectorD cheb2_nodes( double x1, double x2, int num );

VectorD cheb2_weights( int num );

BarFit fit_cheb2_def( Func_Ptr, int n_nodes );

BarFit fit_cheb2( Func_Ptr, int n_nodes, double x1, double x2 );

//plot .c

FILE * gp_open(  );

void gp_term_def( FILE * gp, char * title );

void gp_axes_labels( FILE * gp, char * x_label, char * y_label );

void gp_set_plot( FILE * gp, int num, char ** labels, char ** styles );

void gp_set_logscale( FILE * gp, bool x, bool y );

void gp_prt_carr( FILE * gp, double * xcarr, double * ycarr, int length );

void gp_end( FILE * gp );

typedef struct{
	char * title;
	char ** labels;
	int num;
} PlotInfo;

void plot_2vecs( const VectorD * x, const VectorD * y );

void plot_func( VectorD * domain, Func_Ptr func );

void plot_mult_vecs( const VectorD * x_axis, const VectorD y_values[], 
					int num_vec, PlotInfo );

PlotInfo plot_info_init( int num);

void tmplot_2vecs( const char * cfg_name, const VectorD * xaxis, const VectorD * yaxis);

// roots.c

#define DEF_TOL -1

double root_bis( Par_Func fnc, double x1, double x2, double tol, VectorD * debug);

double root_newt( const Par_Func * fnc, const Par_Func * derivative, 
					double start, double xtol, double ftol, VectorD * debug );


double root_sec( const Par_Func *, double x1, double x2, 
				double xtol, double ftol, VectorD * debug, int * exit_status );

// integration.c

double int_trap( const Par_Func *, double x1, double x2, int n_subint );

double int_simp( const Par_Func *, double x1, double x2, int n_subint );

double legendre(double x, int degree);

double legendre_root( int degree, int root_index, VectorD *debug );

double int_gauss_legendre( Par_Func * fnc, int order, double x1, double x2 );

double cheb_roots( int degree, int root_index );

// tmpplot.c
/*
typedef struct{
	const char *file;
	const char *style;
	const char *label;
	const char *color;
} DatasetDesc;
*/
void tmplot_carrs( const char *cfg_name, const double *xaxis, 
							double **yaxes, int length, int num_yaxes);
/*
DatasetDesc datadesc_init(const char *file, const char *style, const char *label,
	const char *color);
*/


#endif
