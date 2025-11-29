#include <stdio.h>
#include <math.h>
#include "comp_physics.h"


/*
Template: 

double fn( double x, double * unused )
{
	UNUSED(unused);
	return ;
}

Par_Func fn_p = { fn, NULL, 0 };
 */

double f1( double x, double * unused )
{
	UNUSED(unused);
	return x * log(1 + x);
}

Par_Func f1_p = { f1, NULL, 0 };

double f2( double x, double * unused )
{
	UNUSED(unused);
	return x*x*atan(x);
}

Par_Func f2_p = { f2, NULL, 0 };

int main()
{
	
	
	
	exit(EXIT_SUCCESS);
}

