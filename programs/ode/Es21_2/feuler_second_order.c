#include "comp_physics.h"
#include "ode.h"
#include "explot.h"
#include "array.h"
#include "matrix.h"
#include <math.h>

#define NSTEPS 1000

/*
 * The following eq functions are to be understood as equation relating second
 * derivatives to time and function values / derivatives.
 * 
 * Parameters:
 * t         - time, independant variable
 * u[0]      - function value at time t
 * u[1]      - first derivative at time t
 * void *    - optional parameter
 * 
 * Returns:
 * Second derivative
 * 
 */

double eq1(double t, double *u, [[maybe_unused]]void *unused ){
	return -9*u[0]+sin(2*t);}

double eq2(double t, double *u, [[maybe_unused]]void *unused ){
	return 4*u[0]+4*t;}

double eq3(double t, double *u, [[maybe_unused]]void *unused ){
	return -4*u[1]-4*u[0]+t;}
	
double derivative([[maybe_unused]]double t, double *u, [[maybe_unused]]void *unused){
	return u[1];}

/*
 * Following functions are the exact solutions.
 */

double true1(double t){return sin(3*t)/5+2*cos(3*t)+sin(2*t)/5;}

double true2(double t){return 2*cosh(2*t)-t;}

double true3(double t){return (3*t+(double)5/4)*exp(-2*t)+(t-1)/4;}

double trueder1(double t){return 3*cos(3*t)/5-3*2*sin(3*t)+2*cos(2*t)/5;}
		   
double trueder2(double t){return 2*2*sinh(2*t)-1;}
		   
double trueder3(double t){return 3*exp(-2*t)-2*(3*t+(double)5/4)*exp(-2*t)+(double)1/4;}

Array arr_err( const double *x1, const double *x2, int size )
{
	Array arr_of_errs=arr_new(size);
	for(int i=0; i<size; i++)
		ARR(arr_of_errs,i)=x1[i]-x2[i];
	return arr_of_errs;
}

void plot_errors(const Array *t_output, const Array *func_err, 
		const Array *der_err, int counter,  double t0, double f0, double der0, 
	int nsteps, const char *fnc_name)
{
		// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "Error of Euler method for u\'\'(t) = %s, u(%.0f)=%.2lf,"
							" u\'(%.0lf)=%.2lf", 
							fnc_name, t0, f0,t0,der0);
	sprintf(build_name, "euler_2_err%d", counter);
	sprintf(output_name, "euler_2_err%d", counter);

	Global gb_settings = {
		.title=title,
		.xlabel="t",
		.ylabel="Absolute error",
		.logscale=nullptr,
		.build_name=build_name,
		.output_name=output_name
	};
	// Max error over whole vs nsteps plot 
	SeriesSpec err_of_func={
		.x=t_output->data,
		.y=func_err->data,
		.size=nsteps+1,
		.style="l dt 4",
		.label="u_{true}-u",
		.color=nullptr
	};
	// Error over last timestep vs nsteps plot 
	SeriesSpec err_of_der={
		.x=t_output->data,
		.y=der_err->data,
		.size=nsteps+1,
		.style="l dt 4",
		.label="u\'_{true}-u\'",
		.color=nullptr
	};
	SeriesSpec specs[]={err_of_func, err_of_der};
	eplot_multi(specs, 2, &gb_settings);
}

void test_feuler2(Derivative eq, double t0, double tf, double f0, double der0, 
	int nsteps, const char *fnc_name, ScalarFunc true_func, ScalarFunc true_der)
{
	Matrix f_output={0}; /* 	 first row will be for function value, 
		second for its derivative */
	// Static counter to keep track of function calls, for easier file management.
	static int counter=0; counter++;	
	Array t_output={0};
	Derivative diff_equations[]={derivative, eq};
	double init_values[]={f0, der0}; 
	// Call of forward Euler's method
	ode_feuler(t0, tf, init_values, diff_equations, nullptr,
				2, nsteps, &t_output, &f_output);
	// True func values evaluation
	Array true_values=arr_map(&t_output, true_func);
	Array true_der_value=arr_map(&t_output, true_der);
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "Numerical solution of u\'\'(t) = %s, u(%.0f)=%.2lf", 
							fnc_name, t0, f0);
	sprintf(build_name, "euler_2test%d", counter);
	sprintf(output_name, "euler_2test%d", counter);

	Global gb_settings = {
		.title=title,
		.xlabel="t",
		.ylabel="u(t), u\'(t)",
		.logscale=nullptr,
		.build_name=build_name,
		.output_name=output_name
	};
	// Values of numerical solution at timesteps
	SeriesSpec data_desc_f={
		.x=t_output.data,
		.y=f_output.rows[0],
		.size=nsteps+1,
		.style="l",
		.label="u(t) numerical",
		.color=nullptr
	};
	// Values of derivative of numerical solution at timesteps
	SeriesSpec data_desc_der={
		.x=t_output.data,
		.y=f_output.rows[1],
		.size=nsteps+1,
		.style="l",
		.label="u\'(t) numerical",
		.color=nullptr
	};
	// Values of exact solution at timesteps
	SeriesSpec data_desc_true={
		.x=t_output.data,
		.y=true_values.data,
		.size=nsteps+1,
		.style="l dt 4",
		.label="u(t) exact",
		.color=nullptr
	};
	
	SeriesSpec data_desc_true_der={
		.x=t_output.data,
		.y=true_der_value.data,
		.size=nsteps+1,
		.style="l dt 4",
		.label="u\'(t) exact",
		.color=nullptr
	};
	
	SeriesSpec specs[]={data_desc_f, data_desc_der, data_desc_true, data_desc_true_der};
	eplot_multi(specs, 4, &gb_settings);
	
	// Evaluate errors
	Array func_err=arr_err(true_values.data, f_output.rows[0], true_values.size);
	Array der_err=arr_err(true_der_value.data, f_output.rows[1], true_der_value.size);
	
	plot_errors(&t_output, &func_err, &der_err, counter, t0, f0, der0, nsteps,
				fnc_name);
	
	mat_free(&f_output);
	arr_free(&t_output);
	arr_free(&true_values);
	arr_free(&true_der_value);
	arr_free(&func_err);
	arr_free(&der_err);
}

int main()
{
	test_feuler2(eq1, 0, 2*PI, 2, 1, NSTEPS, "sin(2t)-9u", true1, trueder1);
	test_feuler2(eq2, 0, 1.5, 2, -1, NSTEPS, "4y+4t", true2, trueder2);
	test_feuler2(eq3, 0, 4, 1, 0.75, NSTEPS, "t-4y\'-4y", true3, trueder3);
	
	return 0;
}
