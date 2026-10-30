#include "integration.h"
#include <math.h>
#include "rootf.h"
#include "interpolation.h"
#include <limits.h>

#define MAX_LEG 10'000 // maximum degree of a requestable Legendre polynomial
#define ORDER_AUTO -1
#define DOUBLE_EXP_T_TOL DBL_EPSILON*10

double intg_trap( const ParamFuncPtr f, void *params,
			double x1, double x2, int n_subint )
{
	if(n_subint<1)
		raiseErr("n_subint must be at least 1");
	
	const double interval_length = (x2-x1)/n_subint;
	
	double integral = (f(x1,params)+f(x2, params))/2;
	
	for( int i=1; i<n_subint; i++ )
		integral += f(x1 + i*interval_length,params);
	
	integral *= interval_length;
	
	return integral;
}

double intg_simp( const ParamFuncPtr f, void *params,
					double x1, double x2, int n_subint )
{
	const double interval_length = (x2 - x1) / n_subint / 2;
	
	double integral = f(x1,params) + f(x2,params);
	
	double sum_first_points = 0; 
	double sum_second_points = 0;

	for( int i = 1; i <= n_subint; i++ )
		sum_first_points += f(x1 +(2*i - 1)*interval_length, params); 
	
	for( int i = 1; i < n_subint; i++ )
		sum_second_points += f(x1 + 2*i*interval_length,params);
	
	integral += 4 * sum_first_points + 2 * sum_second_points;
	integral *= interval_length / 3;
	
	return integral;	
}


static inline void update_values( double x, int degree, double *values, 
									int stored_degrees )
{
	for( int next_n = stored_degrees; next_n<=degree; next_n++ )
	{
		const double first_addend = ( 2 * next_n - 1 ) * x * values[ next_n - 1 ];
		const double second_addend = ( 1 - next_n ) * values[ next_n - 2 ];
		values[next_n] = (first_addend + second_addend) / next_n;
	}
}

double intg_legendre( double x, int degree )
{
	static double values[MAX_LEG]={1};
	static int stored_degrees = 1;
	static double curr_x = NAN;
	
	if(degree < 0)
	 raiseErr("Degree must be >0\nDegree: %d",degree);
	
	if(degree >= MAX_LEG)
		raiseErr("Requested Legendre polynomial degree must be in [0, %d]\n"
				 "degree = %d",
		MAX_LEG-1,degree);
	
	if( x != curr_x)
	{
		values[1] = x;
		stored_degrees = 2;
		curr_x = x;
	}
		
	if( degree >= stored_degrees )
	{
		update_values( x, degree, values, stored_degrees );
		stored_degrees = degree + 1;
	}

	
	return values[degree];
}

double intg_legendre_wrap(double x, void *params){
	return intg_legendre(x,((int*)params)[0]);}

double intg_legendre_der( double x, int degree )
{
	if(degree == 0)
		return 0;
		
	if(degree<0)
		raiseErr("Degree must be >= 0\n"
		         "degree = %d", degree);
	
	const double num = x*intg_legendre ( x, degree ) - intg_legendre( x, degree - 1 );
	const double den = x*x - 1;
	return degree * num/den;
}

double intg_legendre_der_wrap(double x, void *params){
	return intg_legendre_der(x,((int*)params)[0]);}

static inline double cheb_roots( int degree, int root_index )
{
	const double cos_arg = PI*( 4*root_index - 1 )/( 4*degree + 2);
	return cos(cos_arg);
}

double intg_legendre_root( int degree, int root_index, Array *debug )
{
//	eprint("asked for degree %d and index %d", degree, root_index);
	
	return root_newt( intg_legendre_wrap, &degree, intg_legendre_der_wrap,
						&degree, cheb_roots( degree, root_index), 
				DBL_EPSILON, 1000*DBL_EPSILON, debug );
}

double intg_gauss_legendre_weight( double x, int n )
{
	const double den1 = 1-x*x;
	const double den2 = pown( intg_legendre_der(x, n), 2);
	return 2/den1/den2;
}

double intg_gauss_legendre( ParamFuncPtr f, void *params, int order, 
			double x1, double x2 )
{
	double sum = 0;
	double curr_root = NAN;
	const double half_interval = (x2 - x1)/2;
	const double mid_interval = (x1 + x2)/2;
	for( int i = 1; i <= order; i++)
	{
		curr_root = intg_legendre_root( order, i, NULL );
		const double curr_weight = intg_gauss_legendre_weight( curr_root, order );
		
		const double transformed_root = half_interval * curr_root + mid_interval;
		const double fnc_value = f(transformed_root,params);
		sum += curr_weight * fnc_value;
	}
	return sum * half_interval;
}

static inline double double_exp_real_line_cov(double x){
		return sinh(sinh(x)*PI/2);}

static inline double double_exp_real_line_cov_der(double x){
		return cosh(sinh(x)*PI/2)*cosh(x)*PI/2;}

static inline double real_line_add(double x,ParamFuncPtr f,void *params)
{
	return f(double_exp_real_line_cov(x),params)
					*double_exp_real_line_cov_der(x);}

static inline double double_exp_half_line_cov(double x){
		return exp(PI*sinh(x)/2);}

static inline double double_exp_half_line_cov_der(double x){
		return exp(PI*sinh(x)/2)*PI/2*cosh(x);}

static inline double half_line_add(double x,ParamFuncPtr f,void *params)
{
	return f(double_exp_half_line_cov(x),params)
					*double_exp_half_line_cov_der(x);}

static inline double double_exp_half_line_exp_cov(double x){
		return exp(x-exp(-x));}

static inline double double_exp_half_line_exp_cov_der(double x){
		return exp(x-exp(-x))*(1+exp(-x));}

static inline double half_line_add_exp(double x,ParamFuncPtr f,void *params)
{
	return f(double_exp_half_line_exp_cov(x),params)
					*double_exp_half_line_exp_cov_der(x);}

typedef struct {
    void *params;
	double a;
	double b;
}AddParams;

static inline double double_exp_bounded_cov(double x){
		return tanh(PI/2*sinh(x));}

static inline double double_exp_bounded_cov_der(double x){
		return PI/2*cosh(x)/pow(cosh(PI/2*sinh(x)),2);}

static inline double bounded_add(double t,ParamFuncPtr f,void *params)
{
	AddParams *add_params = (AddParams*)params;

    const double c = (add_params->a+add_params->b)/2;
    const double h = (add_params->b-add_params->a)/2;

    const double u = double_exp_bounded_cov(t);

	return f(c+h*u,add_params->params) * h
					* double_exp_bounded_cov_der(t);
}

static inline double double_exp_findTM(ParamFuncPtr f, void *params,
	double tol, double (*diffeomorphism)(double),
	double (*addend_calc)(double,ParamFuncPtr,void *))
{
	double tM=0;
	double curr_add_plus=0;
	double curr_add_minus=0;
	do{
		tM+=0.1;
		curr_add_plus=fabs(addend_calc(tM,f,params));
		curr_add_minus=fabs(addend_calc(-tM,f,params));
		if(tM>5.7)
			raiseErr("f(x) is not approaching zero fast enough\n"
					 "g(tM) = %g g(-tM) = %g\n"
					 "f(x) = %g f(-x) = %g\n",
						addend_calc(tM,f,params),
						addend_calc(-tM,f,params),
						f(diffeomorphism(tM),params),
						f(diffeomorphism(-tM),params));
		
	}while( curr_add_plus > tol 
		|| curr_add_minus > tol 
		|| isnan(curr_add_plus) || isnan(curr_add_minus));
	return tM;
}

static inline double double_exp_base(ParamFuncPtr f, void *params, int N,
	double tol, double (*diffeomorphism)(double),
	double (*addend_calc)(double,ParamFuncPtr,void *))
{
	const double tM=double_exp_findTM(f, params, tol, diffeomorphism, 
		addend_calc);
	
	double sum_plus=0;
	double sum_minus=0;
	const double stepsize=tM/N;
	
	if(N<=0)
		raiseErr("order N must be greater than 0\n"
				 "order N = %d",N);
	
	for(int i=N; i>0; i--)
	{
		sum_plus+=addend_calc(stepsize*i,f,params);
		sum_minus+=addend_calc(-stepsize*i,f,params);
	}
	return (sum_plus+sum_minus+addend_calc(0,f,params))*stepsize;
}

double intg_double_exp(ParamFuncPtr f, void *params, int N, Dbexp_mode mode,
	double a, double b, double tol)
{
	if(tol==DEF_TOL)
		tol=DOUBLE_EXP_T_TOL;
	else if(tol<=0)
		raiseErr("tol must be more than zero\ntol = %lf",tol);
	
	double (*diffeomorphism)(double);
	double (*addend_calc)(double,ParamFuncPtr,void *);
	AddParams add_params={params,0,0};
	
	switch(mode){
		case REAL_LINE:
			addend_calc=real_line_add;
			diffeomorphism=double_exp_real_line_cov;
			break;
		case HALF_LINE:
			addend_calc=half_line_add;
			diffeomorphism=double_exp_half_line_cov;
			break;
		case HALF_LINE_EXP:
			addend_calc=half_line_add_exp;
			diffeomorphism=double_exp_half_line_exp_cov;
			break;
		case BOUNDED:
			add_params.a=a;
			add_params.b=b;
			params=(void*)&add_params;
			addend_calc=bounded_add;
			diffeomorphism=double_exp_bounded_cov;
	}
	
	return double_exp_base(f, params, N, tol, diffeomorphism, addend_calc);
}
