#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "genutil.h"
#include "arrayf.h"


float averageF( const ArrayF *arr )
{
	float sum = 0;
	for (int i = 0; i < arr->size; i++)
		sum += ARRP(arr,i);
	
	return sum / (float)arr->size;
}

float definition_varianceF( const ArrayF *arr )
{
	if (arr->size == 1)
		return 0;
	
	float ave = averageF( arr );
	float sum = 0;
	for( int i = 0; i < arr->size; i++ )
		sum += ( ARRP(arr,i) - ave )*( ARRP(arr, i) - ave );
		
	return sum / (float)(arr->size - 1);
}

float one_loop_varianceF( const ArrayF *arr )
{
	float sum_squares = 0;
	float sum = 0;
	
	for(int i = 0; i < arr->size; i++)
	{
		sum_squares += ARRP(arr,i)*ARRP(arr,i);
		sum += ARRP(arr,i);
	}
	
	printf("sum = %f\n", sum);
	printf("sum^2/3 = %f\n", - sum*sum / (float)arr->size);
	printf("sum squares = %f\n", sum_squares);
	
	return ( sum_squares - sum*sum / (float)arr->size ) / 
		(float)(arr->size - 1);
}

void print_sample(float values[])
{
	putchar('[');
	for(int i=0; i<3; i++)
	{
		if(i!=0)
			printf("%d+", i);
		printf("%.0e", values[i]);
		if(i!=2)
			fputs(", ", stdout);
	}
	putchar(']');
}

void test_var(float values[])
{
	ArrayF arr=farr_asarr(values,3);
	printf("Of sample ");
	print_sample(values);
	printf( "\nVariance from definition is %.20f\n", definition_varianceF(&arr) );
	printf( "One-loop variance is        %.20f\n\n", one_loop_varianceF(&arr) );
}

int main()
{
	float s1[]={ 1e3f, 1+1e3f, 2+1e3f };
	float s2[]={ 1e6f, 1+1e6f, 2+1e6f };
	float s3[]={ 1e7f, 1+1e7f, 2+1e7f };
	float s4[]={ 1e8f, 1+1e8f, 2+1e8f };
	printf("The sample [1e8f, 1+1e8f, 2+1e8f] gets saved as:\n[");
	for(int i=0; i<3; i++)
	{
		printf("%.0f", s4[i]);
		if(i!=2)
			printf(", ");
	}
	printf("]\n\n");
	
	test_var(s1);
	test_var(s2);
	test_var(s3);
	test_var(s4);
	
	exit(EXIT_SUCCESS);
}
