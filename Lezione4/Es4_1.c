#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "DoubleUtilities.c"

int main()
{
	MatrixDouble testMatrix = readMatD( "startmatrix.txt", 4, 3 );
	printMatDGraph( testMatrix );
	freeMatD( testMatrix );
	
	exit(EXIT_SUCCESS);
}
