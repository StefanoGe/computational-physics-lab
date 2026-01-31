#include "comp_physics.h"
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
	
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "to be added");
	sprintf(build_name, "nonlinear1");
	sprintf(output_name, "nonlinear1");

	Global gb_settings = {
		.title=title,
		.xlabel="x",
		.ylabel="y",
		.logscale=nullptr,
		.build_name=build_name,
		.output_name=output_name
	};
	// Max error over whole vs nsteps plot 
	SeriesSpec curve1={
		.x=sol1.rows[0],
		.y=sol1.rows[1],
		.size=nsteps+1,
		.style="l",
		.label="[x,y](0)=[0.1,0]",
		.color=nullptr
	};
	// Max error over whole vs nsteps plot 
	SeriesSpec curve2={
		.x=sol2.rows[0],
		.y=sol2.rows[1],
		.size=nsteps+1,
		.style="l dt 4",
		.label="[x,y](0)=[0,1.9]",
		.color=nullptr
	};

	SeriesSpec series[]={curve1, curve2};
	eplot_multi(series, 2, &gb_settings);
	
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
	
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "to be added");
	sprintf(build_name, "nonlinear1");
	sprintf(output_name, "nonlinear1");

	Global gb_settings = {
		.title=title,
		.xlabel="x",
		.ylabel="y",
		.logscale=nullptr,
		.build_name=build_name,
		.output_name=output_name
	};

	SeriesSpec curve1={
		.x=sol1.rows[0],
		.y=sol1.rows[1],
		.size=nsteps+1,
		.style="l",
		.label="[x,y](0)=[0.95,0]",
		.color=nullptr
	};

	SeriesSpec curve2={
		.x=sol2.rows[0],
		.y=sol2.rows[1],
		.size=nsteps+1,
		.style="l",
		.label="[x,y](0)=[0,1.05]",
		.color=nullptr
	};
 
	SeriesSpec curve3={
		.x=sol3.rows[0],
		.y=sol3.rows[1],
		.size=nsteps+1,
		.style="l dt 4",
		.label="[x,y](0)=[-2.5,0]",
		.color=nullptr
	};

	SeriesSpec series[]={curve1, curve2, curve3};
	eplot_multi(series, 3, &gb_settings);
	
	Array *arr_to_free[]={&time1, &time2, &time3};
	arr_free_many(arr_to_free, 3);
	Matrix *mat_to_free[]={&sol1, &sol2, &sol3};
	mat_free_many(mat_to_free, 3);
}

int main()
{
	// Problem: how to determine the number of steps more accurately?
	// Scelgo valori sempre più alti, poi confronto le differenze di errori globali
	// - Posso scegliere la "norma" che fa la differenza nel punto finale
	// - Scelgo ogni volta un numero di steps multiplo
	// Reminder: add plot of the circles
	// Idea: add maybe a keep scaling of the axes?
	solve_eq1(1000);
	// Remark: it looks like the two curves gets closer and closer as time 
	// passes without ever touching themselves;
	// indeed a solution is the standard (scaled) parametrization of the unit
	// circle. So we can conclude that the first added is the circling one, 
	// while the second makes the curve converge towards the center.
	solve_eq2(1'000);
	// Remark: we can generalise the previous remark. When the solution is outside
	// the outer circle, the second addend is negative and so the solution gets
	// shrunk onto the outer circle. When it is between the inner and outer circle,
	// it second addend is positive so again the solution converges to the outer circle.
	// When it is inside inner circle, second addend is negative so it makes the solution
	// "collapse" to a point.
	
	return 0;
}
