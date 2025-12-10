#include "comp_physics.h"
#include "ode.h"
#include "explot.h"

double eq1(double t, double *u, [[maybe_unused]]void *unused )
{return -2*t**u;}
/*
void test_feuler(Derivative eq, double t0, double tf, double in_value)
{
	
	
	
	
}
*/

void solve_and_plot(Derivative eq, double t0, double tf, double in_value, int nsteps)
{
	Matrix f_output;
	Array t_output;
	Array initial_value=arr_asarr(&in_value, 1);
	
	ode_feuler( t0, tf, &initial_value, &eq, 1, nsteps, &t_output, &f_output );
	
	
	
}


int main()
{
	
	
	
	
	return 0;
}
