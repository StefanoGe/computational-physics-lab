// vector.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comp_physics.h"

static VectorD NULL_VEC = { NULL, 0, 0 };


VectorD alloc_vecD( int size )
{
	VectorD vec;
	vec.length = 0;
	vec.size = size;
	vec.val = malloc( size * sizeof(double) );
	if (!vec.val)
		raiseErr( "Failed to allocate memory" );
	return vec;
}

VectorD init_vecD(){return alloc_vecD( DEFAULT_VEC_SIZE );}

void appendD(VectorD * vec, double entry)
{
	if (vec->length == vec->size)
	{
		vec->val = realloc( vec->val, vec->size * 2 * sizeof(double) );
		if(!vec->val)
			raiseErr( "Failed to reallocate memory" );
	}
	vec->val[vec->length] = entry;
	vec->length += 1;
}

static inline void free_one_vecD ( VectorD * vec){free( vec->val );}

void _free_all_vecD( VectorD * vec, ... )
{

	va_list list;
	va_start(list, vec);
	
	do{
		free_one_vecD(vec);
		vec = va_arg( list, VectorD * );
		eprint( "%p, %d, %d\n", vec->val, vec->length, vec->size );
	}while( !(vec->val == NULL && vec->length == 0 && vec->size == 0) );
	
	va_end(list);
	
}

VectorD build_vecD ( int length, ... )
{
	va_list list;
	va_start(list, length);
	
	VectorD vec = init_vecD();
	
	for( int i = 0; i < length; i++ )
		appendD( &vec,  va_arg( list, double ) );
	
	va_end(list);
	
	return vec;
}

void print_vecD( const VectorD * vec, char * format, FILE * output )
{
	for( int i = 0; i < vec->length; i ++ )
		fprintf(output,  format, vec->val[i] );
}

void std_print_vecD( const VectorD * vec )
{ 
	putchar('[');
	print_vecD( vec, "%lf, ", stdout );
	putchar(']');
	putchar('\n');
}
