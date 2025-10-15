// Double.c

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "GenUtilities.c"
#include "IntUtilities.c"

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
	
	if (array.val == NULL)
		raiseErr("Errore nell'allocazione della memoria in funzione allocArrD\n");
		
	return array;
}

void freeArrD( ArrayDouble array )
{
	free( array.val );
	array.val = NULL;
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

ArrayDouble linspaceD ( double start, double end, int num_values )
{
	if (num_values < 2)
		raiseErr("You must ask for at least 2 values in function linspaceD");
		
	ArrayDouble array = allocArrD( num_values);
	const int n_increments = num_values - 1;
	for(int i = 0; i <= n_increments; i ++)
		array.val[ i ] = (end - start) * ((double) i / n_increments) + start;
	
	return array;
}

ArrayDouble diffArrD( ArrayDouble array1, ArrayDouble array2)
{
	if (array1.length != array2.length )
		raiseErr("Tried to subtract arrays of different sizes in function "
		"diffArrD\n");
	
	ArrayDouble diffArray = allocArrD( array1.length );
	
	for ( int i = 0; i < array1.length; i++ )
		diffArray.val[i] = array1.val[i] - array2.val[i];
		
	return diffArray;
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


ArrayDouble getRowMatD( MatrixDouble matrix, int row )
{
	ArrayDouble array = allocArrD( matrix.ncols );
	array.val = matrix.val[ row ];
	return array;
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

void cpArrayToRowMatD( MatrixDouble matrix, int row, ArrayDouble array )
{
	ArrayDouble rowMat = getRowMatD( matrix, row );
	cpArrD( array, rowMat );
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

void setValueMatD( MatrixDouble matrix, double value )
{
	for(int i = 0; i < matrix.nrows; i ++)
		for(int j = 0; j < matrix.ncols; j ++)
			matrix.val[i][j] = value;
}

char * firstWS( char * s )
{
	while( (*s != '\0') && ( *s != ' ' ))
		s++;
	return s;
}

ArrayInt parseIndices( const char * select )
{
	int count = 0;
	char * copy = malloc(sizeof(char) * (strlen( select ) + 1) );
	strcpy( copy, select );
	char * token = strtok(copy, " ");
	
	while( token)
	{
		count ++;
		token = strtok(NULL, " ");
	}
	free(copy);

	ArrayInt selIndices = allocArrI( count );
	copy = malloc(sizeof(char) * (strlen( select ) + 1) );
	strcpy( copy, select );
	token = strtok(copy, " ");
	int i = 0;
	while( token )
	{
		printf("%s\n", token);
		selIndices.val[i++] = atoi(token);
		token = strtok(NULL, " ");
	}
	free(copy);
	return selIndices;
}

MatrixDouble sliceMatD( MatrixDouble matrix, char * rowSelect, char * colSelect )
{	
	ArrayInt rowIndices = parseIndices( rowSelect );
	ArrayInt colIndices = parseIndices( colSelect );
	
	printArrIPar(rowIndices);
	printArrIPar(colIndices);
	
	const int newCol = colIndices.length;
	const int newRow = rowIndices.length;
	
	if( newRow > matrix.nrows )
		raiseErr( "Sliced matrix should have no more rows than starting one.\n" );
		
	if( newCol > matrix.ncols )
		raiseErr( "Sliced matrix should have no more cols than starting one.\n" );
	
	MatrixDouble newMatrix = allocMatD( newRow, newCol );
	
	for( int i = 0; i < newRow; i++ )
	{
		for( int j = 0; j < newCol; j++ )
			newMatrix.val[i][j] 
				= matrix.val[ rowIndices.val[i] ][colIndices.val[j]];
	}
	freeArrI( rowIndices );
	freeArrI( colIndices );
	
	return newMatrix;
}

void cpMatToMatD( MatrixDouble source, MatrixDouble dest, 
	int startRow, int startCol )
{
	// add guards later
	
	for( int i = 0; i < source.nrows; i ++ )
		for( int j = 0; j < source.ncols; j++ )
		
			dest.val[startRow + i][startCol + j] = source.val[i][j];
}

void divMat( MatrixDouble matrix, double num, MatrixDouble dest )
{
	// Aggiungi guard
	for( int i = 0; i < matrix.nrows; i ++ )
		for( int j = 0; j < matrix.ncols; j++ )
			dest.val[i][j] = matrix.val[i][j] / num;
}

