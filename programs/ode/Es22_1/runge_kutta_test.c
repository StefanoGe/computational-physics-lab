#include "comp_physics.h"
#include "ode.h"
#include "explot.h"
#include "array.h"
#include "matrix.h"
#include <math.h>

double eq1(double t, double *u, [[maybe_unused]]void *unused )
{return -2*t**u;}

typedef double(*TrueFunc)(double);

double true1(double t)
{
	return 2*exp(-t*t);
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
	const char *fnc_name, TrueFunc fnc)
{
	const int max_nstep=300;
	const int nstep_increment=30;
	const int n_iter=max_nstep/nstep_increment;
	// Variables passed to ode euler solver
	Matrix num_sol={0};
	Array time={0};
	// Global errors
	VectorD errors_ie2=init_vecD();
	VectorD errors_rk4=init_vecD();
	// Number of function evaluation for each iteration of the methods
	VectorD nf_ie2=init_vecD();
	VectorD nf_rk4=init_vecD();

	Array true_values={0};
	VectorD nsteps_vec=init_vecD();
	VectorD one_over_step_size_sq=init_vecD();
	VectorD one_over_step_size_4th=init_vecD();
	
	for(int i=0; i<n_iter; i++)
	{
		// Number of steps
		int nsteps=nstep_increment*(i+1);
		appendD(&nf_ie2, (double)nsteps*2);
		appendD(&nf_rk4, (double)nsteps*4);
		appendD(&one_over_step_size_sq, (double)1/nsteps/nsteps);
		appendD(&one_over_step_size_4th, one_over_step_size_sq.val[i]/nsteps/nsteps);
		appendD(&nsteps_vec, (double)nsteps);
		
		// Call of IE2 method
		ode_ie2( t0, tf, &in_value, &eq, nullptr,
						1, nsteps, &time, &num_sol );
		
		// Computation of true values
		arr_init(&true_values, nsteps+1);
		for(int j=0; j<nsteps+1; j++)
			ARR(true_values,j)=fnc(ARR(time,j));
		
		// Error computation for IE2 method
		appendD(&errors_ie2,  compute_error(true_values.data, num_sol.data,nsteps+1) );
		
		// Call of RK4 method
		ode_rk4( t0, tf, &in_value, &eq, nullptr,
						1, nsteps, &time, &num_sol );
		
		//Error computation for RK4 method
		appendD(&errors_rk4, compute_error(true_values.data, num_sol.data,nsteps+1) );		
		
		arr_free(&time);
		mat_free(&num_sol);
	}
	
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "Convergence of IE2 and RK4 method for u\'(t)=%s", 
							fnc_name);
	sprintf(build_name, "ie2_rk4_conv_study");
	sprintf(output_name, "ie2_rk4_conv_study");

	Global gb_settings = {
		.title=title,
		.xlabel="Number of function evaluations",
		.ylabel="Maximum error",
		.logscale="xy",
		.build_name=build_name,
		.output_name=output_name
	};
	// Max error over whole vs nsteps plot 
	SeriesSpec error_ie2={
		.x=nf_ie2.val,
		.y=errors_ie2.val,
		.size=n_iter,
		.style="lp",
		.label="ie2",
		.color=nullptr
	};
	// Max error over whole vs nsteps plot 
	SeriesSpec error_rk4={
		.x=nf_rk4.val,
		.y=errors_rk4.val,
		.size=n_iter,
		.style="lp",
		.label="rk4",
		.color=nullptr
	};
	// Reference line for second order convergence
	SeriesSpec ref_line2={
		.x=nsteps_vec.val,
		.y=one_over_step_size_sq.val,
		.size=n_iter,
		.style="l dt 4",
		.label="2^{nd} order convergence",
		.color=nullptr
	};
	// Reference line for fourth order convergence
	SeriesSpec ref_line4={
		.x=nsteps_vec.val,
		.y=one_over_step_size_4th.val,
		.size=n_iter,
		.style="l dt 4",
		.label="4^{th} order convergence",
		.color=nullptr
	};

	SeriesSpec series[]={error_ie2, error_rk4, ref_line2, ref_line4};
	eplot_multi(series, 4, &gb_settings);
			
	free_vecD(&errors_ie2, &errors_rk4, &nf_ie2, &nf_rk4, &nsteps_vec,
				&one_over_step_size_4th, &one_over_step_size_sq);
}

int main()
{
	error_study(eq1, 0, 2, 2, "-2tu", true1);
	
	return 0;
}
