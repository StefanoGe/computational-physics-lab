#include "integration.h"
#include <math.h>
#include "rootf.h"
#include "interpolation.h"
#include <limits.h>

#define MAX_LEG 10'000
#define ORDER_AUTO -1
#define DOUBLE_EXP_T_TOL DBL_EPSILON*5

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
					*double_exp_real_line_cov_der(x);
}



double intg_double_exp_real_line(ParamFuncPtr f, void *params, int N)
{
	double (diffeomorphism)(double,ParamFuncPtr,void*);
	
	double tM=0.1;
	while( fabs(real_line_add(tM,f,params)) > DOUBLE_EXP_T_TOL 
		|| fabs(real_line_add(-tM,f,params)) > DOUBLE_EXP_T_TOL )
	{
		if(tM>5.7)
			raiseErr("f(x) is not approaching zero fast enough\n"
					 "g(tM) = %g g(-tM) = %g\n"
					 "f(x) = %g f(-x) = %g\n",
						real_line_add(tM,f,params),
						real_line_add(-tM,f,params),
						f(double_exp_real_line_cov(tM),params),
						f(double_exp_real_line_cov(-tM),params));
/*	
		eprint(	"%lf\n"
				"g(tM) = %g g(-tM) = %g\n",
				 tM,
				 real_line_add(tM,f,params),
					real_line_add(-tM,f,params));
*/
		tM+=0.1;
	}
	
//	eprint("%lf\n",tM);
	
	double sum=0;
	const double stepsize=tM/N;
	
	if(N<=0)
		raiseErr("order N must be greater than 0\n"
				 "order N = %d",N);
	
	for(int i=-N; i<=N; i++)
	{
//		eprint("%lf, %lf\n",double_exp_real_line_cov(i*stepsize),
//						double_exp_real_line_cov_der(i*stepsize));
		sum+=real_line_add(stepsize*i,f,params);
/*
		eprint( "sum = %lf addend = %lf\n"
				"x = %lf f(x) = %lf",sum,real_line_add(stepsize*i,f,params),
					double_exp_real_line_cov(stepsize*i),
					f(double_exp_real_line_cov(stepsize*i),params));
*/
	}
	
	return sum*stepsize;
}
