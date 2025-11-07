// vector.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comp_physics.h"

const VectorD NULL_VEC = { NULL, 0, 0 };

static inline int necessary_size(int length)
{
	int size = DEFAULT_VEC_SIZE;
	for( ; size < length; size *=2 );
	return size;
}

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

VectorD init_vec_length( int length )
{
	VectorD vec = alloc_vecD( necessary_size(length) );
	vec.length = length;
	return vec;
}

void appendD(VectorD * vec, double entry)
{
	if (vec->length == vec->size)
	{
		vec->val = realloc( vec->val, vec->size * 2 * sizeof(double) );
		if(!vec->val)
			raiseErr( "Failed to reallocate memory" );
		vec->size *= 2;
	}
	vec->val[vec->length] = entry;
	++vec->length;
}

static inline void free_one_vecD ( VectorD * vec){free( vec->val );}

static inline bool are_same_vecD( const VectorD * vec1, const VectorD * vec2 )
{
	return vec1->val == vec2->val && vec1->length == vec2->length 
		&& vec1->size == vec2->size;
}

void _free_all_vecD( VectorD * vec, ... )
{
	va_list list;
	va_start(list, vec);
	
	do{
		free_one_vecD(vec);
		vec = va_arg( list, VectorD * );
		//eprint( "%p, %d, %d\n", (void *)vec->val, vec->length, vec->size );
	}while( !(are_same_vecD( vec, &NULL_VEC )) );
	
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

VectorD vec_cp( const VectorD * source )
{
	VectorD new_vec = alloc_vecD( source -> size );
	new_vec.length = source ->length;
	for( int i = 0; i < new_vec.length; i++ )
		new_vec.val[i] = source -> val [i];
	return new_vec;
}

VectorD vec_range( double x1, double x2, int num )
{
	if( x2 < x1 )
		raiseErr( "given x2 < x1" );
	const double spacing = (x2- x1) / (num - 1);
	VectorD vec = init_vecD();
	for( int i = 0; i < num; i++ )
		appendD( &vec, i * spacing + x1 );
	
	return vec;
}
