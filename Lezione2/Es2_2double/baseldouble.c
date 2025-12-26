#include <stdlib.h>
#include <math.h>
#include "genutil.h"
#include "array.h"
#include "explot.h"

#define SAMPLE_STEP 1.2
#define SAMPLE_SIZE 110
#define INT_OVERF_GUARD (long int)pow(2,40)

double sumSquaresDouble( long int maxIndex )
{
	double sum = 0;
	for ( int i = 1; i <= maxIndex; i ++ )
		sum += 1.0 / (double)i /(double) i;
		
	return sum;
}

double revSumSquaresDouble( long int maxIndex )
{
	double sum = 0;
	for ( long int i = maxIndex; i > 0; i -- )
		sum += 1.0 / (double)i /(double) i;
		
	return sum;
}


void sumSquaresAnalysis( double sample_step, int sample_size )
{
	// Setup the numbers of iterations for the plotting sample
	const long int max_int = INT_OVERF_GUARD;
	long int *num_of_iterations=malloc(sizeof(long int)*sample_size);
	for(int i=0; i<sample_size; i++)
	{	
		long int iter=(long int)pow(sample_step,i);
		if( iter>=max_int )
			raiseErr("Risk of integer overflow %ld at iter %d", iter, i);
		num_of_iterations[i]=iter;
	}
	
	Array xaxis = arr_new(sample_size);
	for(int i=0; i<sample_size; i++)
		ARR(xaxis,i)=(double)num_of_iterations[i];
	
	Array err_incr = arr_new(sample_size);
	Array err_decr = arr_new(sample_size);
	Array ref_line = arr_new(sample_size);
	
	const double true_value = (double)(PI*PI) / 6 ;
	
	for (int i = 0; i < sample_size; i++)
	{
		ARR(err_incr,i) = true_value - sumSquaresDouble( num_of_iterations[i] );
		ARR(err_decr,i) = true_value - revSumSquaresDouble( num_of_iterations[i] );
		ARR(ref_line,i) = 1 / (double)num_of_iterations[i];
		eprint("done %4d/110 %ld", i, num_of_iterations[i]);
	}

	Global gb_settings = {
		.title="Approximation error of Basel problem - 64bit precision",
		.xlabel="N",
		.ylabel="Error",
		.logscale="xy",
		.build_name="baseldouble",
		.output_name="basel1_2_1double",
		.key="top right",
		.additional_set=nullptr
	};

	SeriesSpec increasing={
		.x=xaxis.data,
		.y=err_incr.data,
		.size=sample_size,
		.style="l",
		.label="Sum with decreasing addends",
		.color=nullptr
	};
	
	SeriesSpec decreasing={
		.x=xaxis.data,
		.y=err_decr.data,
		.size=sample_size,
		.style="l",
		.label="Sum with increasing addends",
		.color=nullptr
	};
	
	SeriesSpec reference={
		.x=xaxis.data,
		.y=ref_line.data,
		.size=sample_size,
		.style="l dt 2",
		.label="1/x reference line",
		.color=nullptr
	};
	
	SeriesSpec series[]={increasing, decreasing, reference};
	eplot_multi(series, 3, &gb_settings);

	Array *arr_to_free[]={&err_incr, &err_decr, &ref_line, &xaxis};
	arr_free_many(arr_to_free, 4);
	free(num_of_iterations);
}


int main ()
{
	printf("Exact value is %lf\n", PI*PI/6);
	sumSquaresAnalysis(SAMPLE_STEP, SAMPLE_SIZE);
	exit(EXIT_SUCCESS);
}
