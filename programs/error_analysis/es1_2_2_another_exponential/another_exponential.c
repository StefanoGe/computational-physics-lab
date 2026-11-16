#include <stdio.h>
#include "genutil.h"
#include <math.h>

double f1(double x)
{
	return (exp(x)-1)/x;
}

double f1_maclaurin(double x, int n)
{
	double sum=1;
	for(int i=0; i<n; i++)
	{
		sum *= x/(n+1-i);
		sum += 1;
	}
	return sum;
}

int main()
{
	for(int i=3; i<=16; i++)
	{
		double x = pow(10,-i);
		printf("%.20lf %.20lf\n", f1(x), f1_maclaurin(x, 4));
	}
	
	return 0;
}
