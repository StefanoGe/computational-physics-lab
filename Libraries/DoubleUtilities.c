// Double.c

#ifndef SG_DOUBLE_UTILITIES_C
#define SG_DOUBLE_UTILITIES_C

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "GenUtilities.c"
#include "IntUtilities.c"
#include <stdbool.h>

#define CREATE_MAT NULL
#define NULL_ARR (ArrayDouble){NULL,0}
#define NULL_MAT (MatrixDouble){NULL, 0,0}

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

void printArrDPar (ArrayDouble array, char * format)
{
	putchar('[');
	for(int i = 0; i < array.length; i ++)
		printf(format, array.val[i]);
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
	FILE * file = openFile( filename, "r" );
	
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

void cpArrayToColMatD( ArrayDouble array, MatrixDouble dest, int col,
		bool freeArr )
{
	if( dest.nrows != array.length )
		raiseErr( "In function cpArrayToColMatD array and dest must be"
				  "of same length!");
	
	for(int i = 0; i < array.length; i++)
		dest.val[i][col] = array.val[i];
	
	if(freeArr)
		freeArrD( array );
}

MatrixDouble transposeMatD( MatrixDouble matrix)
{
	MatrixDouble returnMatrix = allocMatD( matrix.ncols, matrix.nrows );

	for(int i = 0; i < matrix.nrows; i++)
		for(int j = 0; j < matrix.ncols; j++)
			returnMatrix.val[j][i] = matrix.val[i][j];
			
	return returnMatrix;
}

void printDatMatD( MatrixDouble matrix, char * filename, char * format, 
	bool doTranspose)
{
	FILE * fileDat = openFile( filename, "w" );
	
	if(doTranspose)
	{
		for(int i = 0; i < matrix.ncols; i++)
		{
			for(int j = 0; j < matrix.nrows; j++)
				fprintf( fileDat, format, matrix.val[j][i]);
			
			putc('\n', fileDat);
		}
	}else{
		for(int i = 0; i < matrix.nrows; i++)
		{
			for(int j = 0; j < matrix.ncols; j++)
				fprintf( fileDat, format, matrix.val[i][j]);
			
			putc('\n', fileDat);
		}
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


void diagMatD( MatrixDouble matrix, double value )
{
	if( matrix.nrows != matrix.ncols )
		raiseErr("matrix must be square! In func diagMatD\n");
		
	int const dim = matrix.nrows;
	
	for( int i = 0; i < dim; i++ )
		for( int j = 0; j < dim; j++ )
			matrix.val[ i ][ j ] = ( i == j ) ? value : 0;
}

MatrixDouble matMultD( MatrixDouble mat1, MatrixDouble mat2, MatrixDouble * dest)
{	
	MatrixDouble tempMat = allocMatD( mat1.nrows, mat2.ncols );
	setValueMatD( tempMat, 0 );
	
	for( int i = 0; i < mat1.nrows; i++ )
		for( int j = 0; j < mat2.ncols; j++ )
			for( int p = 0; p < mat1.ncols; p++ )
			
				tempMat.val[i][j] += mat1.val[i][p]* mat2.val[p][j];
	/*
	if(bool)
	{
		freeMatD(mat1);
		freeMatD(mat2);
	}
	*/
	
	// Input CREATE_MAT as dest in order to use as dest a newly allocated matrix.
	if( dest == CREATE_MAT )
		return tempMat;
	
	freeMatD( *dest );
	dest -> val = tempMat.val;
	dest -> ncols = tempMat.ncols;
	dest -> nrows = tempMat.nrows;
	
	return *dest;
}

MatrixDouble asColMatrix( ArrayDouble arr, bool destroySource )
{
	MatrixDouble matrix = allocMatD( arr.length, 1 );
	
	for( int i = 0; i < arr.length; i++ )
		matrix.val[i][0] = arr.val[i];
		
	if(destroySource)
		freeArrD( arr );
		
	return matrix;
}

void freeAllArrD( ArrayDouble arr, ... )
{
	va_list list;
	va_start(list, arr);
	
	do{
		freeArrD(arr);
		arr = va_arg( list, ArrayDouble );
		//fprintf(stderr, "a\n");
	}while( !(arr.val == NULL && arr.length == 0) );
	
	va_end(list);
}

void freeAllMatD( MatrixDouble matrix, ... )
{
	va_list list;
	va_start(list, matrix);
	
	do{
		freeMatD(matrix);
		matrix = va_arg( list, MatrixDouble );
		//fprintf(stderr, "a\n");
	}while( !(matrix.val == NULL && matrix.ncols == 0) );
	
	va_end(list);
}

void exchange_rows( MatrixDouble matrix, int row1, int row2 )
{
	if( row1 >= matrix.nrows || row2 >= matrix.nrows )
		raiseErr( "In func exchange_rows:\n"
					"row indeces must be less than matrix.nrows."
					"row1 = %d, row2 = %d", row1, row2 );
	
	if( row1 < 0 || row2 < 0 )
		raiseErr( "In func exchange_rows:\n"
					"indeces should be more than 0: "
					"row1 = %d, row2 = %d", row1, row2 );
	
	if( row1 == row2 )
		return;
	
	double temp = 0;
	const int ncols = matrix.ncols;
	
	for( int col = 0; col < ncols; col++ )
	{
		temp = matrix.val[row1][col];
		matrix.val[row1][col] = matrix.val[row2][col];
		matrix.val[row2][col] = temp;
	}
}

ArrayDouble mat_vec_mult( MatrixDouble matrix, ArrayDouble array, 
ArrayDouble * dest )
{
	if( matrix.ncols != array.length )
		raiseErr( "In func mat_vec_mult:"
					"matrix.mcols != array.length\n"
					"matrix.ncols = %d, matrix.length = %d\n", 
					matrix.ncols, array.length);
					
	MatrixDouble col_array = asColMatrix( array, false );

	matMultD(matrix, col_array, &col_array);
	
	ArrayDouble final_arr = allocArrD( array.length );
	
	for( int i = 0; i < array.length; i ++)
		final_arr.val[i] = col_array.val[i][0];
	
	freeMatD( col_array );
	
	if( dest == NULL )
		return final_arr;
	
	freeArrD( *dest );
	dest -> val = final_arr.val;
	dest ->length = final_arr.length;
	
	return *dest;
	
	
}

MatrixDouble mat_diffD( MatrixDouble mat1, MatrixDouble mat2, MatrixDouble * dest)
{	
	MatrixDouble tempMat = allocMatD( mat1.nrows, mat2.ncols );
	setValueMatD( tempMat, 0 );
	
	for( int i = 0; i < mat1.nrows; i++ )
		for( int j = 0; j < mat2.ncols; j++ )
			
				tempMat.val[i][j] = mat1.val[i][j] - mat2.val[i][j];
	
	// Input CREATE_MAT as dest in order to use as dest a newly allocated matrix.
	if( dest == CREATE_MAT )
		return tempMat;
	
	freeMatD( *dest );
	dest -> val = tempMat.val;
	dest -> ncols = tempMat.ncols;
	dest -> nrows = tempMat.nrows;
	
	return *dest;
}



#endif
