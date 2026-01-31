#include <stdlib.h>
#include <math.h>
#include "genutil.h"
#include "arrayf.h"
#include "array.h"
#include "explot.h"

#define MAX_INDEX 10000

float sumSquaresFloat( int maxIndex )
{
	float sum = 0;
	for ( int i = 1; i <= maxIndex; i ++ )
		sum += 1.f / (float) (i * i);
		
	return sum;
}

float revSumSquaresFloat( int maxIndex )
{
	float sum = 0;
	for ( int i = maxIndex; i > 0; i -- )
		sum += 1.f / (float) (i * i);
		
	return sum;
}


void sumSquaresAnalysis( int maxIndex )
{
	ArrayF err_incr = arrf_new(maxIndex+1);
	ArrayF err_decr = arrf_new(maxIndex+1);
	const float true_value = (float)(PI*PI) / 6 ;
	
	ARR(err_incr,0)=0;
	ARR(err_decr,0)=0;
	
	for (int i = 1; i <= maxIndex; i++)
	{
		ARR(err_incr,i) = true_value - sumSquaresFloat( i );
		ARR(err_decr,i) = true_value - revSumSquaresFloat( i ) ;
	}
	
	Array err_incrD = arr_from_float(&err_incr);
	Array err_decrD = arr_from_float(&err_decr);
	
	Array xaxis = arr_new(maxIndex+1);
	for(int i=0; i<=maxIndex; i++)
		ARR(xaxis,i)=(double)i;

	Global gb_settings = {
		.title="Approximation error of Basel problem - 32bit precision",
		.xlabel="N",
		.ylabel="Error",
		.logscale=nullptr,
		.build_name="baselfloat",
		.output_name="basel1_2_1float",
		.key="top right",
		.additional_set="set xrange [2000:10000]\n"
	};

	SeriesSpec increasing={
		.x=xaxis.data,
		.y=err_incrD.data,
		.size=maxIndex+1,
		.style="l",
		.label="Sum with decreasing addends",
		.color=nullptr
	};
	
	SeriesSpec decreasing={
		.x=xaxis.data,
		.y=err_decrD.data,
		.size=maxIndex+1,
		.style="l",
		.label="Sum with increasing addends",
		.color=nullptr
	};
	
	SeriesSpec series[]={increasing, decreasing};
	eplot_multi(series, 2, &gb_settings);

}


int main ()
{
	printf("Exact value is %lf\n", PI*PI/6);
	sumSquaresAnalysis(MAX_INDEX);
	exit(EXIT_SUCCESS);
}

void testSumSquare()
{
	printf("Ordine normale:\n");
	
	printf("%f\n", sumSquaresFloat( 30 ) );
	
	printf("Ordine reversed:\n");
	
	printf("%f\n", revSumSquaresFloat( 30 ) );
	
	printf( "Valore vero:\n%f", PI*PI / 6 );
}
