#include "comp_physics.h"
#include "ode.h"
#include "explot.h"

double eq1(double t, double *u, [[maybe_unused]]void *unused )
{return -2*t**u;}

double e2(double t, double *u, [[maybe_unused]]void *unused )
{return *u+t;}

double eq3(double t, double *u, [[maybe_unused]]void *unused )
{return t*t/(*u)/(1+t*t*t);}



void solve_and_plot(Derivative eq, double t0, double tf, double in_value, 
		int nsteps, const char *fnc_name)
{
	Matrix f_output;
	Array t_output;
	Array initial_value=arr_asarr(&in_value, 1);
	static int counter=0;
	counter++;
	
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "Numerical solution of u\'(t) = %s, u(%.0f)=%.0lf", 
							fnc_name, t0, in_value);
	sprintf(build_name, "euler_test%d", counter);
	sprintf(output_name, "euler_test%d", counter);

	Global gb_settings = {
		.title=title,
		.xlabel="t",
		.ylabel="u(t)",
		.logscale=nullptr,
		.build_name=build_name,
		.output_name=output_name
	};
	
	
	
	ode_feuler( t0, tf, &initial_value, &eq, 1, nsteps, &t_output, &f_output );

	eplot_2carr(t_output.data, f_output.data, nsteps+1, &gb_settings,
			&NULL_DATADESC);
}

void test_feuler(Derivative eq, double t0, double tf, double in_value, 
	const char *fnc_name)
{
	solve_and_plot(eq, t0, tf, in_value, 320, fnc_name);
	
	
	
}


int main()
{
	
	test_feuler(eq1, 0, 2, 2, "-2tu");
	
	
	return 0;
}
