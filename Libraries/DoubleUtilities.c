// Double.c

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "GenUtilities.c"
#include "string.h"

typedef struct _arrayDouble
{
	double * val;
	int length;
	
} ArrayDouble;

typedef struct _MatrixDouble
{
	double ** val;
	int nrows;
	int ncols;
} MatrixDouble;

ArrayDouble allocArrD ( int length)
{
	ArrayDouble array;
	array.length = length;
	array.val = (double *)malloc(length * sizeof(double));
	if (array.val == NULL) {
		printf("Errore nell'allocazione della memoria in funzione allocArrD\n");
		exit(1);
    }
	return array;
}

MatrixDouble allocMatD ( int nrows, int ncols )
{
	MatrixDouble matrix;
	matrix.nrows = nrows;
	matrix.ncols = ncols;
	
	matrix.val = (double **)malloc( nrows * sizeof(double*) );
	for( int i = 0; i < nrows; i++ )
		matrix.val[i] = (double *)malloc(ncols *sizeof(double));
		
	return matrix;
}

void printMatDGraph( MatrixDouble matrix )
{
	putchar('[');
	for( int i = 0; i < matrix.nrows; i++ )
	{
		if( i != 0 )
			putchar('|');
		
		for( int j = 0; j < matrix.ncols - 1; j++ )
			printf("% 10lf,", matrix.val[i][j]);
		
		printf("% 10lf", matrix.val[i][matrix.ncols - 1]);
		if( i + 1 < matrix.nrows )
			printf("|\n");	
	}
	putchar(']');
	putchar('\n');
}

MatrixDouble readMatD(char * filename, int nrows, int ncols)
{
	FILE * file = fopen( filename, "r" );
	
	MatrixDouble matrix = allocMatD( nrows, ncols );
	for( int i = 0; i < nrows; i ++ )
		for( int j = 0; j < ncols; j++ )
			fscanf( file, "%lf", matrix.val[i] + j );
	
	fclose(file);
	return matrix;
}

void freeMatD( MatrixDouble matrix )
{
	for(int i = 0; i < matrix.nrows; i++)
	{
		free(matrix.val[i]);
		matrix.val[i] = NULL;
	}
	free(matrix.val);
	matrix.val = NULL;
}

ArrayDouble buildArrD ( int length, ... )
{
	va_list list;
	va_start(list, length);
	
	ArrayDouble array = allocArrD( length );
	
	for( int i = 0; i < length; i++ )
		array.val[ i ] = va_arg( list, double );
		
	va_end(list);
	
	return array;
}

void freeArrD( ArrayDouble array )
{
	free( array.val );
	array.val = NULL;
}

void printArrDPar (ArrayDouble array)
{
	putchar('[');
	for(int i = 0; i < array.length; i ++)
		printf("%.10lf, ", array.val[i]);
	printf("]\n");
}

void setValueArrD(ArrayDouble array, double value)
{
	for( int i = 0; i < array.length; i++ )
		array.val[i] = value;
}


ArrayDouble diffArrD( ArrayDouble array1, ArrayDouble array2)
{
	if (array1.length != array2.length )
	{
		printf("Tried to subtract arrays of different sizes in function "
		"diffArrD\n");
		exit(EXIT_FAILURE);
	}
	
	ArrayDouble diffArray = allocArrD( array1.length );
	
	for ( int i = 0; i < array1.length; i++ )
		diffArray.val[i] = array1.val[i] - array2.val[i];
		
	return diffArray;
}

ArrayDouble getRowMatD( MatrixDouble matrix, int row )
{
	ArrayDouble array = allocArrD( matrix.ncols );
	array.val = matrix.val[ row ];
	return array;
}

void cpArrD( ArrayDouble source, ArrayDouble destination )
{
	if( source.length != destination.length )
	{
		printf( "In function cpArrD source and destination must be of"
		"the same length! [%d; %d]\n", source.length, destination.length );
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < source.length; i ++)
		destination.val[i] = source.val[i];
}

void passArrayToRowMatD( MatrixDouble matrix, int row, ArrayDouble array )
{
	if (matrix.ncols != array.length)
	{
		printf("In function passArrayToRowMatD matrix and array"
		"must be of same height! [%d, %d]\n", matrix.ncols, array.length);
		exit(EXIT_FAILURE);
	}
	
	matrix.val[row] = array.val;
}

MatrixDouble transposeMatD( MatrixDouble matrix)
{
	MatrixDouble returnMatrix = allocMatD( matrix.ncols, matrix.nrows );

	for(int i = 0; i < matrix.nrows; i++)
		for(int j = 0; j < matrix.ncols; j++)
			returnMatrix.val[j][i] = matrix.val[i][j];
			
	return returnMatrix;
}

void printDatMatD( MatrixDouble matrix, char * filename, char * format )
{
	FILE * fileDat = openFile( filename, "w" );
	
	for(int i = 0; i < matrix.nrows; i++)
	{
		for(int j = 0; j < matrix.ncols; j++)
			fprintf( fileDat, format, matrix.val[i][j]);
		
		putc('\n', fileDat);
	}
	fclose(fileDat);
}
