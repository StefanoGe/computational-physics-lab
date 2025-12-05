// Es2_2.c

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# define A 111
# define B 1130
# define C 3000
# define X_0 11.0f/2
# define X_1 61.0f/11  
# define REQUIRED_INDEX 34


float numInSequence( int n )
{
	if( n == 0 )
	{
		return X_0;
	}
	
	if( n == 1 )
	{
		return X_1;
	}
	
	float first_div = C / numInSequence( n - 2 );
	
	printf( "First_div of required index %d is: %f\n", n, first_div );
	
	float first_subt = B - first_div;
	
	printf( "First_sub of required index %d is: %f\n", n, first_subt );
	
	float second_div = first_subt/ numInSequence( n - 1 );

	printf( "Second_div of required index %d is: %f\n", n, second_div );

	return A - second_div;
}

float seqCalc( float minus1, float minus2, int n )
{
	float first_div = C / minus2;
	
	printf( "First_div of required index %d is: %f\n", n, first_div );
	
	float first_subt = B - first_div;
	
	printf( "First_sub of required index %d is: %f\n", n, first_subt );
	
	float second_div = first_subt/ minus1;

	printf( "Second_div of required index %d is: %f\n", n, second_div );

	float returnValue = A - second_div;
	
	printf("X of required index %d is: %f\n", n, returnValue);

	return returnValue;
}

void sequence()
{
	float seq[REQUIRED_INDEX + 1] = {0};
	seq[0] = X_0;
	seq[1] = X_1;
	
	for (int i = 2; i <= REQUIRED_INDEX; i ++)
		seq[i] = seqCalc( seq[i-1], seq[ i - 2], i );
}

int main()
{

	sequence();
		
	exit(EXIT_SUCCESS);
}
/*
for(int i = 0; i < REQUIRED_INDEX; i ++ )
{
	printf("X of required index %d is: %f\n", 
		i, numInSequence(i));
}
*/
