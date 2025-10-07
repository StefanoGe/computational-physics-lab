#include <iostream>
#include <cmath>
#include <vector>
#include <utility>

double approxExpDank( double x, int truncationIndex)
{
	double sum{1};
	
	for( int i = truncationIndex; i > 0; i-- )
	{
		sum *= x /i;
		sum += 1;
	}
	
	return sum;
}

double pow ( double x, int n )
{
	
	return n == 0 ? 1 : x * pow( x, n-1 );
}

int fact( int n)
{
	if ( n < 0 )
		return 1;
	
	return n == 0 ? 1 : n * fact(n-1);
}


double approxExpBasic ( double x, int truncIndex )
{
	if ( truncIndex == 0 )
		return 1;
		
	return approxExpBasic( x, truncIndex - 1 ) + pow( x, truncIndex ) / fact(truncIndex);
}



void testDouble( double x, int nDigits = 10 )
{
	std::cout.precision(nDigits);
	std::cout << "Il valore ottenuto e`: " << x << "\n";
}

double absError(double x, int truncIndex)
{
	
	return fabs( approxExpDank( x, truncIndex ) - exp(x) );
}

void testExps( double testValue )
{	
	std::cout << "Basic:\n";
	
	testDouble ( approxExpBasic( testValue, 10 ) );
	
	std::cout << "Dank:\n";
	
	testDouble ( approxExpDank( testValue, 10 ) );
	
	std::cout << "Abs err dank\n";
	
	testDouble ( absError( testValue, 10 ) );
}

void printVec ( std::vector<double> vec )
{
	for( double const& elem : vec  )
		std::cout << elem << "\n";
}

void printPairVec ( std::vector<std::pair<double, double>> vec )
{
	for( auto const& pair : vec  )
		std::cout << "[ " << pair.first << " , " << pair.second << " ]\n";
}

void absErrNGrowth ( double x, int lastN )
{

	std::vector<double> absErrVec{};

	for ( int i = 0; i <= lastN; i++ )
		absErrVec.emplace_back( absError( x, i ) );
		
	std::vector<double> relErrVec{};
	
	for ( int i = 0; i <= lastN; i++ )
		relErrVec.emplace_back( absErrVec.data()[i] * fact( i + 1 ) / pow( x, i + 1 ) );

	printVec( absErrVec );
	printVec( relErrVec );

}

void absErrXGrowth ( int numXValues, int N )
{
	std::vector<double> xvalues;

	for( int i = 0; i <= numXValues; i++ )
		xvalues.emplace_back( static_cast<double>( i ) / numXValues );
	
	std::vector<std::pair<double, double>> absErrVec{};

	for ( const auto& x : xvalues )
		absErrVec.emplace_back( x, absError( x, N ) );
	
	std::vector<std::pair<double, double>> relErrVec{};
	
	for ( const auto& pair : absErrVec )
		relErrVec.emplace_back( pair.first, pair.second * fact( N + 1 ) / pow( pair.first, N + 1 ) );

	printPairVec( absErrVec );
	printPairVec( relErrVec );

}

int main()
{	
	
	absErrXGrowth( 20, 2 );
	
	return 0;
}

void testHello()
{	std::cout << "Hello world!" << std::endl; }



