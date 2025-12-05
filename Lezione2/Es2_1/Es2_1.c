// Esercizio 1 della lezione 2

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define MAX_INDEX 10000

FILE * openFile( char * fileName, char * mode )
{
	FILE * file;
	file = fopen(fileName, mode);
	if (file == NULL) {
		printf("The file %s is not opened.", fileName);
    }
    return file;
}

typedef struct _arrayFloat
{
	float * val;
	int length;
	
} ArrayFloat;

ArrayFloat allocArrF ( int length)
{
	ArrayFloat array;
	array.length = length;
	array.val = (float *)malloc(length * sizeof(float));
	return array;
}

void freeArrF( ArrayFloat array )
{
	free( array.val );
	array.val = NULL;
}

float sumSquaresFloat( int maxIndex )
{
	float sum = 0;
	for ( int i = 1; i <= maxIndex; i ++ )
		sum += 1.f / (float) (i * i);
		
	return sum;
}

float revSumSquaresFloat( int maxIndex )
{
	float sum = 0;
	for ( int i = maxIndex; i > 0; i -- )
		sum += 1.f / (float) (i * i);
		
	return sum;
}

void printArrFPar (ArrayFloat array)
{
	putchar('[');
	for(int i = 0; i < array.length; i ++)
		printf("%.10f, ", array.val[i]);
	printf("]\n");
}

void setValueArrF(ArrayFloat array, float value)
{
	for( int i = 0; i < array.length; i++ )
		array.val[i] = value;
}


ArrayFloat diffArrF( ArrayFloat array1, ArrayFloat array2)
{
	if (array1.length != array2.length )
	{
		printf("Tried to subtract arrays of different sizes in function "
		"diffArrF\n");
		exit(EXIT_FAILURE);
	}
	
	ArrayFloat diffArray = allocArrF( array1.length );
	
	for ( int i = 0; i < array1.length; i++ )
		diffArray.val[i] = array1.val[i] - array2.val[i];
		
	return diffArray;
}

void sumSquaresAnalysis( int maxIndex )
{
	ArrayFloat partSums = allocArrF(maxIndex + 1);
	ArrayFloat revPartSums = allocArrF(maxIndex + 1);
	
	partSums.val[0] = 0;
	revPartSums.val[0] = 0;
	
	for (int i = 1; i <= maxIndex; i++)
	{
		partSums.val[i] = partSums.val[ i -  1 ] += 1/(float) (i * i);
		revPartSums.val[i] = revSumSquaresFloat( i ) ;
	}
	
	ArrayFloat trueArr = allocArrF(maxIndex + 1);
	setValueArrF ( trueArr, (float)(PI*PI) / 6 );
	
	ArrayFloat truncErr = diffArrF( trueArr, partSums );
	ArrayFloat revTruncErr = diffArrF( trueArr, revPartSums );
	
	//printArrFPar( truncErr );
	//printArrFPar( revTruncErr );
	
	FILE * output = openFile( "data.dat", "w" );
	
	for(int i = 1; i <= maxIndex; i++)
		fprintf( output,"%d %f %f\n",
			i, truncErr.val[i], revTruncErr.val[i]);
	
	freeArrF( partSums );
	freeArrF( revPartSums );
	freeArrF( trueArr );
	freeArrF( truncErr ) ;
	freeArrF( revTruncErr );
	fclose( output );
	
	printf("Done!\n");
}


int main ()
{
	sumSquaresAnalysis(MAX_INDEX);
	exit(EXIT_SUCCESS);
}

void testArray()
{
	ArrayFloat testArray = allocArrF( 10 );
	for (int i = 0; i < testArray.length; i++)
		testArray.val[i] = (float) (i*i);
	
	printArrFPar( testArray );
}

void testSumSquare()
{
	printf("Ordine normale:\n");
	
	printf("%f\n", sumSquaresFloat( 30 ) );
	
	printf("Ordine reversed:\n");
	
	printf("%f\n", revSumSquaresFloat( 30 ) );
	
	printf( "Valore vero:\n%f", PI*PI / 6 );
}
