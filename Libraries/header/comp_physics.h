// My general header

# ifndef SG_COMP_PHYSICS_HEADER_H
# define SG_COMP_PHYSICS_HEADER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>

#define PI 3.1415926535897932384626433

#define CREATE_MAT NULL
#define NULL_ARR (ArrayDouble){NULL,0}
#define NULL_MAT (MatrixDouble){NULL, 0,0}
#define NULL_INFO (PlotInfo){NULL, NULL,0}

#define UNUSED(expr) do { (void)(expr); } while(0)

#define raiseErr( msg, ...) \
	do { \
		fprintf(stderr, "[%s: line %d in %s]: " msg, __FILE__, __LINE__, __func__ \
		__VA_OPT__(,) __VA_ARGS__); \
		exit(EXIT_FAILURE); \
	} while (0)

#define eprint(a, ...) fprintf( stderr, a "\n" __VA_OPT__(,) __VA_ARGS__ )

typedef double(*Func_Ptr)(double);

typedef double(*ParamFuncPtr) (double, double *);

typedef struct {
	ParamFuncPtr param_func_ptr;
	double * params;
	int nparams;
} Par_Func;

double evaluate( const Par_Func * self, double x );

FILE * openFile( char * fileName, char * mode );

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

typedef struct MatrixDouble MatrixDouble;

// double * __m_get( MatrixDouble * mat, int nrow, int ncol );

struct MatrixDouble
{
	double ** val;
	int nrows;
	int ncols;
	//double * (*get)( MatrixDouble * self, int nrow, int ncol );
};

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

void gp_set_plot( FILE * gp, int num, char ** labels );

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


#endif
