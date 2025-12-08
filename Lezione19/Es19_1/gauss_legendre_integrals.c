#include "comp_physics.h"
#include <math.h>

#define MAX_N_NODES 40

double f1( double x, void *unused ){
	UNUSED(unused);
	return exp(-4*x);}
	
double f2( double x, void *unused ){
	UNUSED(unused);
	return exp(-9*x*x);}
	
double f3( double x, void *unused ){
	UNUSED(unused);
	return 1/cosh(x);}

double f4( double x, void *unused ){
	UNUSED(unused);
	return 1/(1 + 9*x*x);}
	
double f5( double x, void *unused ){
	UNUSED(unused);
	return x*x*sin(8*x);}

void computation_result_analysis( int max_n_nodes, ParamFuncPtr func, double true_value,
	double x1, double x2)
{
	Par_Func func_p = { func, 0, 0 };
	VectorD evaluated_integrals = init_vecD();
	VectorD true_value_deviation = init_vecD();
	int n_nodes[ 1000 ];
	int n_estimations = max_n_nodes / 4;
	
	
	for( int i = 0; i < n_estimations; i++ )
	{
		n_nodes[i] = ( i+1 ) * 4;
		appendD(&evaluated_integrals, int_gauss_legendre( &func_p, n_nodes[i],
								x1, x2));
		appendD(&true_value_deviation, fabs(true_value - evaluated_integrals.val[i]) );
	}
	
	std_print_vecD(&evaluated_integrals);
	
	VectorD n_nodesD = init_vecD();
	for(int i = 0; i < n_estimations; i++)
		appendD( &n_nodesD, (double) n_nodes[i]);
	
	tmplot_2vecs( "gauss_legendre", &n_nodesD, &true_value_deviation );
	// Ideally, the plot should show all the functions at the same time.
	
	free_vecD( &evaluated_integrals, &true_value_deviation, &n_nodesD );
}


int main()
{
	computation_result_analysis( MAX_N_NODES, f1, sinh(4)/2, -1, 1 );
	computation_result_analysis( MAX_N_NODES, f2, sqrt(PI) * erf(3)/3, -1, 1 );
	computation_result_analysis( MAX_N_NODES, f3, 2*atan( sinh(1) ), -1, 1 );
	computation_result_analysis( MAX_N_NODES, f4, 2*atan(3)/3, -1 ,1 );
	computation_result_analysis( MAX_N_NODES, f5, -3*PI*PI/32, PI/2 , PI ); 
	// Here it got the exact value in the second to last calculation
	
	// Additional work idea: calculate and verify numerically
	
	return 0;
}
