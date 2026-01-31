// GenericUtilities.c
#include "genutil.h"
#define MAX_FACTORIAL 12

FILE * openFile( const char * fileName, const char * mode )
{
	FILE * file;
	file = fopen(fileName, mode);
	if (file == NULL) 
		raiseErr( "File %s cannot be opened\n", fileName );
		
    return file;
}

void mypause ( void ) 
{ 
  printf ( "Press [Enter] to continue . . ." );
  fflush ( stdout );
  getchar();
} 

/*
 * factorial function
 * 
 * memorises data with static variables, so recalling the function for
 * a value already computed does NOT compute it again.
 * 
 * To prevent integer overflow, maximum value that can be calculated is 
 * factorial(MAX_FACTORIAL-1)
 * 
 */

int factorial( int n )
{
	static int stored_values[MAX_FACTORIAL] = {1,1};
	static int count = 2; // stored_values[0] and [1] are already valid

    if(n < 0)
		raiseErr("negative argument %d", n);
	
	if(n >= MAX_FACTORIAL)
		raiseErr("maximum that can be evaluated is set to %d", MAX_FACTORIAL);
	
	if (n >= count){
		for(int i = count; i<=n; i++ )
			stored_values[i] = stored_values[i-1]*i;
		count=n+1;
	}
	
	return stored_values[n];
}

double pown( double x, int n )
{
	double product = 1.0;
	for( int i = 0; i < n; i++ )
		product *=x;
	
	return product;
}

void carr_print_inline(const double *arr, int size, const char *format, 
	bool parentheses, bool linebreak)
{
	if (parentheses)
		fputs("[ ", stdout);
	
	for(int i = 0; i < size; i++ )
	{
		if(i != 0)
			fputs( ", ", stdout);

		printf(format, arr[i]);
		
	}
	if (parentheses)
		fputs( " ]", stdout);
	if(linebreak)
		putchar('\n');
}





