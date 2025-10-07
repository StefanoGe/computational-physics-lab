// Float.c

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
