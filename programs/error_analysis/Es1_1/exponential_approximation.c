#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "genutil.h"
#include "matrix.h"
#include "array.h"
#include "explot.h"

#define MAX_N 4
#define N_INTERVALS 100

double naiveExp( double x, int N )
{
	double sum = 0;
	
	for(int i=0; i<=N; i++)
		sum += pown(x, i) / factorial(i);
		
	return sum;
}

double smartExp( double x, int N )
{
	double sum = 1;
	for( int i = 0; i < N; i++ )
	{
		sum *= x/( N - i );
		sum += 1;
	}
	return sum;
}

double absoluteError( double x, int N )
{
	return fabs( smartExp( x, N ) - exp( x ) );
}

void analysis ( int max_N, int n_intervals )
{
	const int n_xvalues = n_intervals + 1;
	Array xAxis = arr_linspace( 0, 1, n_xvalues);
	Matrix errs = mat_new(max_N, n_xvalues);
	Matrix ref = mat_new(max_N, n_xvalues);

	for( int i = 1; i <= max_N; i++ )
	{
		for(int j = 1; j < n_xvalues; j++)
		{
			MAT(errs, i-1,j)=absoluteError( ARR(xAxis,j), i);
//			printf("x = %lf - exp = %lf - err = %lf\n", 
//				ARR(xAxis,j), smartExp(ARR(xAxis,j), i), MAT(errs,i-1,j));
			MAT(ref,i-1,j-1)=pown(ARR(xAxis,j-1),i+1)/factorial(i+1);
			
			//eprint("%lf %d",MAT(ref,i-1,j-1), factorial(i));
		}
	}
		
	Global gb_settings = {
		.title="Approximation error of exponential",
		.xlabel="x",
		.ylabel="Truncation error",
		.logscale="xy",
		.build_name="approx_exp",
		.output_name="approx_exp1_0",
		.key="bottom right font \",12\""
	};

	SeriesSpec n1={
		.x=xAxis.data,
		.y=errs.rows[0],
		.size=n_xvalues,
		.style="l",
		.label="N=1",
		.color=nullptr
	};
	
	SeriesSpec n2={
		.x=xAxis.data,
		.y=errs.rows[1],
		.size=n_xvalues,
		.style="l",
		.label="N=2",
		.color=nullptr
	};
	
	SeriesSpec n3={
		.x=xAxis.data,
		.y=errs.rows[2],
		.size=n_xvalues,
		.style="l",
		.label="N=3",
		.color=nullptr
	};
	
	SeriesSpec n4={
		.x=xAxis.data,
		.y=errs.rows[3],
		.size=n_xvalues,
		.style="l",
		.label="N=4",
		.color=nullptr
	};

	SeriesSpec ref1={
		.x=xAxis.data,
		.y=ref.rows[0],
		.size=n_xvalues,
		.style="l dt 4",
		.label="x^{2} / 2!",
		.color=nullptr
	};
	
	SeriesSpec ref2={
		.x=xAxis.data,
		.y=ref.rows[1],
		.size=n_xvalues,
		.style="l dt 4",
		.label="x^{3} / 3!",
		.color=nullptr
	};
	
	SeriesSpec ref3={
		.x=xAxis.data,
		.y=ref.rows[2],
		.size=n_xvalues,
		.style="l dt 4",
		.label="x^{4} / 4!",
		.color=nullptr
	};
	
	SeriesSpec ref4={
		.x=xAxis.data,
		.y=ref.rows[3],
		.size=n_xvalues,
		.style="l dt 4",
		.label="x^{5} / 5!",
		.color=nullptr
	};

	SeriesSpec series[]={n1,ref1,n2,ref2,n3,ref3,n4,ref4};
	eplot_multi(series, 8, &gb_settings);


	arr_free(&xAxis);
	mat_free(&errs);
}


int main()
{
	analysis(MAX_N, N_INTERVALS);
	
	return 0;
}

/* 
 * Remarks:
 * the difference between e^x and the approximating series
 * is e^x - eapprox, which is asymptotic to x^n+1/n+1! as x approaches 0.
 * So the graphs look almost the same in this limit, but start to diverge when
 * when approaching 1.
 */

