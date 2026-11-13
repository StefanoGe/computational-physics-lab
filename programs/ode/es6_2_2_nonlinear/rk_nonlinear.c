#include "ode.h"
#include "explot.h"
#include "array.h"
#include "matrix.h"
#include <math.h>

const double t0=0;
const double tf=10;

double eqx1([[maybe_unused]] double t, double *u, [[maybe_unused]]void *unused )
{
	const double x=u[0]; const double y=u[1];
	return -4*y+x*(1-x*x-y*y);
}

double eqy1([[maybe_unused]] double t, double *u, [[maybe_unused]]void *unused )
{
	const double x=u[0]; const double y=u[1];
	return 4*x+y*(1-x*x-y*y);
}

double eqx2([[maybe_unused]] double t, double *u, [[maybe_unused]]void *unused )
{
	const double x=u[0]; const double y=u[1];
	return -4*y-x/4*(1-x*x-y*y)*(4-x*x-y*y);
}

double eqy2([[maybe_unused]] double t, double *u, [[maybe_unused]]void *unused )
{
	const double x=u[0]; const double y=u[1];
	return 4*x-y/4*(1-x*x-y*y)*(4-x*x-y*y);
}

typedef struct {
	double maxdiff;
	double lastdiff;
} DiffWithExact;

DiffWithExact compare_solution_max_distance(const Matrix *sol1, const Matrix *sol2, 
	int ratio)
{
	DiffWithExact difference = {0};
	for(int i=0; i<sol1->ncols; i++)
	{
		difference.maxdiff=fmax(difference.maxdiff, 
			fabs(MATP(sol1,0,i)-MATP(sol2,0,i*ratio)));
		difference.maxdiff=fmax(difference.maxdiff, 
			fabs(MATP(sol1,1,i)-MATP(sol2,1,i*ratio)));
	}
	difference.lastdiff=fabs(MATP(sol1,0,sol1->ncols-1)-MATP(sol2,0,sol2->ncols-1));
	difference.lastdiff=fmax(difference.lastdiff, 
		fabs(MATP(sol1,1,sol1->ncols-1)-MATP(sol2,1,sol2->ncols-1)));
	return difference;
}

void analyze_the_accuracy(const Matrix *sol,
		const double in_value[], const Derivative eqs[], int nsteps)
{
	Array time_check={0};
	Matrix sol_check={0};
	ode_rk4( t0, tf, in_value, eqs, nullptr, 2, nsteps*2, &time_check,
		&sol_check );
	
	DiffWithExact difference = compare_solution_max_distance(sol,
		&sol_check, 2);
	
	printf(ANSI_GREEN"Difference between solution in %d steps and solution in \n"
			"%d steps for curve starting in [%g, %g]: \n"
			ANSI_YELLOW"    At last time step: "ANSI_RED"%g\n"
			ANSI_YELLOW"    Maximum error: "ANSI_RED"%g\n"ANSI_RESET"\n",
			nsteps, nsteps*2, in_value[0], in_value[1],
			difference.lastdiff, difference.maxdiff);
		
	arr_free(&time_check);
	mat_free(&sol_check);
}

void solve_eq1(int nsteps)
{
	const double in_value1[]={0.1,0};
	const double in_value2[]={0,1.9};
	Matrix sol1={0};
	Array time1={0};

	Matrix sol2={0};
	Array time2={0};
	const Derivative eqs[]={eqx1, eqy1};
	
	ode_rk4( t0, tf, in_value1, eqs, nullptr, 2, nsteps, &time1, &sol1 );
	ode_rk4( t0, tf, in_value2, eqs, nullptr, 2, nsteps, &time2, &sol2 );
	
	// Retake with double points
	
	analyze_the_accuracy(&sol1, in_value1, eqs, nsteps);
	analyze_the_accuracy(&sol2, in_value2, eqs, nsteps);
	
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "to be added");
	sprintf(build_name, "nonlinear1");
	sprintf(output_name, "nonlinear1");

	Global gb_settings = eplot_init_global(title, "x", "y", nullptr, build_name,
		output_name, "", NULL);
	// Max error over whole vs nsteps plot 
	SeriesSpec curve1=eplot_init_series(sol1.rows[0],sol1.rows[1],nsteps+1,
		"l","[x,y](0)=[0.1,0]",nullptr);
	// Max error over whole vs nsteps plot 
	SeriesSpec curve2=eplot_init_series(sol2.rows[0],sol2.rows[1],nsteps+1,
		"l dt 4","[x,y](0)=[0,1.9]",nullptr);

	SeriesSpec series[]={curve1, curve2};
	eplot_multi_options(series, 2, &gb_settings, false);
	
	puts("\n------------------------------------------------\n");
	
	Array *arr_to_free[]={&time1, &time2};
	arr_free_many(arr_to_free, 2);
	Matrix *mat_to_free[]={&sol1, &sol2};
	mat_free_many(mat_to_free, 2);
}

void solve_eq2(int nsteps)
{
	const double in_value1[]={0.95,0};
	const double in_value2[]={0,1.05};
	const double in_value3[]={-2.5,0};
	Matrix sol1={0};
	Array time1={0};
	Matrix sol2={0};
	Array time2={0};
	Matrix sol3={0};
	Array time3={0};
	const Derivative eqs[]={eqx2, eqy2};
	
	ode_rk4( t0, tf, in_value1, eqs, nullptr, 2, nsteps, &time1, &sol1 );
	ode_rk4( t0, tf, in_value2, eqs, nullptr, 2, nsteps, &time2, &sol2 );
	ode_rk4( t0, tf, in_value3, eqs, nullptr, 2, nsteps, &time3, &sol3 );
	
	analyze_the_accuracy(&sol1, in_value1, eqs, nsteps);
	analyze_the_accuracy(&sol2, in_value2, eqs, nsteps);
	analyze_the_accuracy(&sol3, in_value3, eqs, nsteps);
	
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "to be added");
	sprintf(build_name, "nonlinear1");
	sprintf(output_name, "nonlinear1");

	Global gb_settings = eplot_init_global(title,"x","y",nullptr,build_name,
		output_name, "", NULL);

	SeriesSpec curve1=eplot_init_series(sol1.rows[0],sol1.rows[1],nsteps+1,"l",
		"[x,y](0)=[0.95,0]",nullptr);

	SeriesSpec curve2=eplot_init_series(sol2.rows[0], sol2.rows[1], nsteps+1,
		"l", "[x,y](0)=[0,1.05]",nullptr);
 
	SeriesSpec curve3=eplot_init_series(sol3.rows[0],sol3.rows[1],nsteps+1,"l dt 4",
		"[x,y](0)=[-2.5,0]",nullptr);

	SeriesSpec series[]={curve1, curve2, curve3};
	eplot_multi_options(series, 3, &gb_settings, false);
	
	Array *arr_to_free[]={&time1, &time2, &time3};
	arr_free_many(arr_to_free, 3);
	Matrix *mat_to_free[]={&sol1, &sol2, &sol3};
	mat_free_many(mat_to_free, 3);
}

int main()
{
	solve_eq1(1000);
	// Remark: it looks like the two curves gets closer and closer as time 
	// passes without ever touching themselves;
	// indeed a solution is the standard (scaled) parametrization of the unit
	// circle. So we can conclude that the first added is the circling one, 
	// while the second makes the curve converge towards the center.
	solve_eq2(1000);
	// Remark: we can generalise the previous remark. When the solution is outside
	// the outer circle, the second addend is negative and so the solution gets
	// shrunk onto the outer circle. When it is between the inner and outer circle,
	// it second addend is positive so again the solution converges to the outer circle.
	// When it is inside inner circle, second addend is negative so it makes the solution
	// "collapse" to a point.
	
	return 0;
}
