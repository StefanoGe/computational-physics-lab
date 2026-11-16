// Es 3_1double.c
#include "array.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


double average( const Array *arr )
{
	double sum = 0;
	for (int i = 0; i < arr->size; i++)
		sum += ARRP(arr,i);
	
	return sum / (double) arr->size;
}

double definition_variance( const Array *arr )
{
	if (arr->size == 1)
		return 0;
	
	double ave = average( arr );
	double sum = 0;
	for( int i = 0; i < arr->size; i++ )
		sum += ( ARRP(arr,i) - ave )*( ARRP(arr,i) - ave );
		
	return sum / (arr->size - 1);
}

double one_loop_variance( const Array *arr )
{
	double sum_squares = 0;
	double sum = 0;
	
	for (int i = 0; i < arr->size; i++)
	{
		sum_squares += ARRP(arr,i)*ARRP(arr,i);
		sum += ARRP(arr,i);
	}
	
//	printf("sum = %f\n", sum);
//	printf("sum squares = %f\n", sum_squares);
	
	return ( sum_squares - sum*sum / arr->size ) / 
		(arr->size - 1);
}

void test_var(double values[])
{
	Array arr=arr_asarr(values,3);
	
	printf("Of sample [");
	
	for(int i=0; i<3; i++)
	{
		if(i!=0)
			printf("%d+", i);
		printf("%.0e", ARR(arr, i));
		if(i!=2)
			fputs(", ", stdout);
	}
	puts("]");
	printf( "Variance from definition is %.20lf\n", definition_variance(&arr) );
	printf( "One-loop variance is        %.20lf\n\n", one_loop_variance(&arr) );
}


int main()
{
	double sample1[] = {1e3, 1+1e3, 2+1e3};
	double sample2[] = {1e6, 1+1e6, 2+1e6};
	double sample3[] = {1e7, 1+1e7, 2+1e7};
	double sample4[] = {1e8, 1+1e8, 2+1e8};
	
	test_var(sample1);
	test_var(sample2);
	test_var(sample3);
	test_var(sample4);
	
	exit(EXIT_SUCCESS);
}


