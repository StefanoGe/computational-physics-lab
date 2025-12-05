#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <comp_physics.h>

double e_sin(double x)
{
	return exp(sin(x-1));
}

double polyn_approx( double x, ArrayDouble parameters )
{
	double sum = 0;
	if ( x < 1 )
		for( int n = parameters.length - 1; n >= 0; n-- )
			sum += parameters.val[n]*pown(x, n);
	else
		for( int n = 0; n < parameters.length; n++ )
			sum += parameters.val[n]*pown(x, n);
	
	return sum;
}

double test_polyn_approx( double x )
{
	ArrayDouble parameters = buildArrD( 2, 0.5, -1.0 );
	
	return polyn_approx( x, parameters );
}

void plot_to_file( double x_min, double x_max, int count, char * filename, 
						Func_Ptr func)
{
	ArrayDouble x_data = linspaceD( x_min, x_max, count );
	ArrayDouble y_data = build_from( x_data, func );	

	FILE * file = openFile( filename, "w" );
	for(int i = 0; i < count; i++)
		fprintf(file,"%.10lf %.10lf\n", x_data.val[i], y_data.val[i] );
	fclose(file);
}

void approx_e_sin()
{
	ArrayDouble x_data = linspaceD( 2 * PI / 60, 1, 60 );
	ArrayDouble y_data = build_from( x_data, e_sin );
	
	ArrayDouble estimated_parameters = linear_ls_pol_fitting( 
														x_data, y_data, 6 );
	
	ArrayDouble x_data_plot = linspaceD(0, 1, 100 );
	MatrixDouble plot_data = allocMatD(100, 2);
	for( int i = 0; i < plot_data.nrows; i++ )
	{
		plot_data.val[i][0] = x_data_plot.val[i];
		plot_data.val[i][1] = polyn_approx( x_data_plot.val[i], estimated_parameters );
	}
	
	printDatMatD(plot_data, "data/expsin_fitted_plot.dat" ,"%lf ",false);
	
	printArrDPar( estimated_parameters, "%lf, " );
	freeAllArrD( x_data, y_data, estimated_parameters, x_data_plot, NULL_ARR);
	freeMatD(plot_data);
}

int main()
{
	approx_e_sin();
	
	//plot_to_file( 0, 1, 100, "data/test_data.dat", test_polyn_approx );
	
	return 0;
}





