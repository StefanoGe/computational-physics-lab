#include "comp_physics.h"
#include "ode.h"
#include "explot.h"
#include "array.h"
#include "matrix.h"
#include <math.h>

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
	return -4*y+x/4*(1-x*x-y*y)*(4-x*x-y*y);
}

double eqy2([[maybe_unused]] double t, double *u, [[maybe_unused]]void *unused )
{
	const double x=u[0]; const double y=u[1];
	return 4*x+y/4*(1-x*x-y*y)*(4-x*x-y*y);
}

void solve_eq1(Derivative eqx, Derivative eqy, double t0, double tf, double x0,
			double y0, int nsteps)
{
	Matrix sol={0};
	Array time={0};
	const Derivative eqs[]={eqx, eqy};
	const double in_value[]={x0,y0};
	
	static int counter=0;
	counter++;
	
	ode_rk4( t0, tf, in_value, eqs, nullptr, 2, nsteps, &time, &sol );
	
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "to be added");
	sprintf(build_name, "nonlinear%d", counter);
	sprintf(output_name, "nonlinear%d", counter);

	Global gb_settings = {
		.title=title,
		.xlabel="x",
		.ylabel="y",
		.logscale=nullptr,
		.build_name=build_name,
		.output_name=output_name
	};
	// Max error over whole vs nsteps plot 
	SeriesSpec curve={
		.x=sol.rows[0],
		.y=sol.rows[1],
		.size=nsteps+1,
		.style="l",
		.label="boh",
		.color=nullptr
	};

	SeriesSpec series[]={curve};
	eplot_multi(series, 1, &gb_settings);
	
}

void solve_eq2(Derivative eqx, Derivative eqy, double t0, double tf, double x0,
			double y0, int nsteps)
{
	Matrix sol={0};
	Array time={0};
	const Derivative eqs[]={eqx, eqy};
	const double in_value[]={x0,y0};
	
	static int counter=0;
	counter++;
	
	ode_rk4( t0, tf, in_value, eqs, nullptr, 2, nsteps, &time, &sol );
	
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "to be added");
	sprintf(build_name, "nonlinear%d", counter);
	sprintf(output_name, "nonlinear%d", counter);

	Global gb_settings = {
		.title=title,
		.xlabel="x",
		.ylabel="y",
		.logscale=nullptr,
		.build_name=build_name,
		.output_name=output_name
	};
	// Max error over whole vs nsteps plot 
	SeriesSpec curve={
		.x=sol.rows[0],
		.y=sol.rows[1],
		.size=nsteps+1,
		.style="l",
		.label="boh",
		.color=nullptr
	};

	SeriesSpec series[]={curve};
	eplot_multi(series, 1, &gb_settings);
	
}

int main()
{
	solve_eq(eqx1, eqy1, 0, 10, 0.1, 0, 1000);
	solve_eq(eqx2, eqy2, 0, 10, 0.95, 0, 1000);
	
	return 0;
}
