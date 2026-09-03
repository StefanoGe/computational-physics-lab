#include "genutil.h"
#include "ode.h"
#include "explot.h"
#include "array.h"
#include "matrix.h"
#include <math.h>

#define N_FUNCS 3

double eq1(double t, double *u, [[maybe_unused]]void *unused )
{return -2*t**u;}

double eq2(double t, double *u, [[maybe_unused]]void *unused )
{return *u+t;}

double eq3(double t, double *u, [[maybe_unused]]void *unused )
{return t*t/(*u)/(1+t*t*t);}

int powint(int base, int exp)
{
	int prod=1;
	for(int i=0; i<exp; i++)
		prod*=base;
	return prod;
}

double true1(double t)
{
	return 2*exp(-t*t);
}

double true2(double t)
{
	return -1-t+3*exp(t);
}

double true3(double t)
{
	return sqrt(1 + 2*log(1+t*t*t)/3);
}
/*
void solve()
{
	// Call of Euler's method.
	ode_feuler( t0, tf, &in_value, &eq, NULL, 1, nsteps, &t_output, &f_output );
}
*/
void plot(
		int nsteps, char fnc_names[][100], FuncPtr funcs[],
		Matrix *f_values, Array *t_values)
{
	// Static counter to keep track of function calls, for easier file management.
	static int counter=0; counter++;
	Array true_values[N_FUNCS]={0};
	Array normalized_times[N_FUNCS]={0};
	for(int i=0; i<N_FUNCS; i++)
	{
		arr_cp(t_values+i,normalized_times+i);
		arr_scale(normalized_times+i,1.0/ARRP(t_values+i,nsteps));
	}
	eprint("ok2");
	for(int i=0; i<N_FUNCS; i++)
		true_values[i]=arr_map(t_values+i,funcs[i]);
	eprint("ok3");
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	char labels[N_FUNCS*2][100];
	sprintf(title, "ODEs solutions from Euler's method");
	sprintf(build_name, "euler_test%d", counter);
	sprintf(output_name, "euler_test%d", counter);
	SeriesSpec *series;
	SAFE_ALLOC(series,N_FUNCS*2);

	Global gb_settings = eplot_init_global(title,"t","u(t)",NULL,build_name,
	output_name,"top right outside",NULL);
	for(int i=0; i<N_FUNCS; i++)
	{
		sprintf(labels[i*2],"u' = %s - exact solution",fnc_names[i]);
		series[i*2]=eplot_init_series(normalized_times[i].data,f_values[i].data,
			nsteps+1,"l",labels[2*i],NULL);
		sprintf(labels[i*2+1],"-    numerical solution");
		series[i*2+1]=eplot_init_series(normalized_times[i].data,true_values[i].data,
			nsteps+1,"l dt 3",labels[2*i+1],NULL);
	}
	/*
	// Values of exact solution at timesteps
	SeriesSpec data_desc_true={
		.x=t_output.data,
		.y=true_values.data,
		.size=nsteps+1,
		.style="lines",
		.label="Exact function",
		.color=nullptr
	};
	*/
	/*
	// Values of numerical solution at timesteps
	SeriesSpec data_desc_est={
		.x=t_output.data,
		.y=f_output.data,
		.size=nsteps+1,
		.style="points",
		.label="Numerical solution",
		.color=nullptr
	};
	*/
	
	eplot_multi(series, 2*N_FUNCS, &gb_settings);
	
	for(int i=0; i<3; i++){
		arr_free(true_values+i);
		arr_free(normalized_times+i);
	}
}

double compute_error(double *true_arr, double *estimate, int length)
{
	double max=0;
	double temp=0;
	for(int i=0; i<length; i++)
	{
		temp=fabs(true_arr[i]-estimate[i]);
		if(temp>=max)
			max=temp;
//		eprint("True value: %.6lf ; Estimate: %.6lf; max = %lf", true_arr[i], estimate[i], max);
	}
//	mypause();
	return max;
}

void error_study(Derivative eq, double t0, double tf, double in_value,
	const char *fnc_name, int max_k, FuncPtr fnc)
{
	// Variables passed to ode euler solver
	Matrix f_output={0};
	Array t_output={0};
	// Variables for error analysis
	Array errors={0};
	Array errorsf={0};
	Array nsteps_vec={0};
	Array approx_stepsize={0};
	Array true_values={0};
	// Static counter to keep track of function calls, for easier file management.
	static int counter=0;
	counter++;
	
	for(int i=0; i<=max_k-2; i++)
	{
		// Number of steps per Euler's method call
		int nsteps=10*powint(2,i+2);
		arr_append(&nsteps_vec, (double)nsteps);
		arr_append(&approx_stepsize, 1/(double)nsteps);
		
		// Call of Euler's method
		ode_feuler( t0, tf, &in_value, &eq, nullptr,
						1, nsteps, &t_output, &f_output );
		
		// Error computation
		arr_init(&true_values, nsteps+1);
		for(int j=0; j<nsteps+1; j++)
			ARR(true_values,j)=fnc(t_output.data[j]);
		arr_append(&errors,  compute_error(true_values.data, f_output.data,nsteps+1) );
		arr_append(&errorsf, fabs( ARR(f_output,nsteps) - true_values.data[nsteps]) );
		
		arr_free(&t_output);
		mat_free(&f_output);
	}
	
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "Convergence of Euler method for u(t)=%s, u(%.0lf)=%.0lf", 
							fnc_name, t0, in_value);
	sprintf(build_name, "euler_conv_study%d", counter);
	sprintf(output_name, "euler_conv_study%d", counter);

	Global gb_settings =eplot_init_global(title,"Number of steps","Absolute error",
		"xy",build_name,output_name,"","");
	// Max error over whole vs nsteps plot 
	SeriesSpec error_overall=eplot_init_series(nsteps_vec.data,errors.data,
		max_k-1,"lp pt 2","||u-v||_{" GP_INFINITY "}",NULL);
	// Error over last timestep vs nsteps plot 
	SeriesSpec error_final={nsteps_vec.data,errorsf.data,max_k-1,"lp pt 1",
		"Error at last time timestep",NULL};
	// Reference line for first order convergence
	SeriesSpec ref_line={nsteps_vec.data,approx_stepsize.data,max_k-1,"lp",
		"First order convergence reference line",NULL};
	
	SeriesSpec series[]={error_overall, error_final, ref_line};
	eplot_multi(series, 3, &gb_settings);
	
	Array *arr_to_free[]={&errors,&errorsf,&nsteps_vec,&approx_stepsize,
		&true_values};
	arr_free_many(arr_to_free,5);
}

/*
void test_feuler(Derivative eq, double t0, double tf, double in_value, 
	const char *fnc_name, FuncPtr true_func)
{
	solve_and_plot(eq, t0, tf, in_value, 320, fnc_name, true_func);
	error_study(eq, t0, tf, in_value, fnc_name, 10, true_func);
}
*/

int main()
{
	Matrix *f_values;
	SAFE_ALLOC(f_values,3);
	Array *t_values;
	SAFE_ALLOC(t_values,3);
	
	const double init_values[]={2,2,1};
	Derivative eqs[]={eq1,eq2,eq3};
	FuncPtr true_funcs[]={true1,true2,true3};
	
	ode_feuler(0,2,init_values  ,eqs  ,NULL,1,320,t_values,  f_values);
	ode_feuler(0,1,init_values+1,eqs+1,NULL,1,320,t_values+1,f_values+1);
	ode_feuler(0,3,init_values+2,eqs+2,NULL,1,320,t_values+2,f_values+2);
	eprint("ok");
	char fnames[][100]={"-2tu","u+t","t^2/u(1+t^3)"};
	plot(320,fnames,true_funcs,f_values,t_values);
/*
	test_feuler(eq1, 0, 2, 2, "", true1);
	test_feuler(eq2, 0, 1, 2, "", true2); 
	// Remark: inf error is the same as error on final timestep. The error
	// then grows steadily.
	test_feuler(eq3, 0, 3, 1, "", true3);
*/
	for(int i=0; i<3; i++)
	{
		SAFE_FREE(f_values);
		SAFE_FREE(t_values);
	}
	
	return 0;
}
