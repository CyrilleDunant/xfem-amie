//
// C++ Interface: random generator
//
// Description: 
//
//
// Author: Alain Giorla
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "random.h"


namespace Mu
{

RandomNumber::RandomNumber() 
{
// 	reset() ;
}

double RandomNumber::normal() 
{
	double y = exponential(1.) ;
	double floor = std::exp(-(y-1.)*(y-1.)/2) ;
	while(uniform() > floor) 
	{
			y = exponential(1.) ;
			floor = std::exp(-(y-1.)*(y-1.)/2) ;
	}
	double sgn = uniform() ;
	if(sgn > 0.5)
		return -y ;
	return y ;
}

double RandomNumber::logistic(double a, double b)
{
	double u = uniform() ;
	return a + b * std::log(u / (1-u)) ;
}

double RandomNumber::xhi2(int a)
{
	double xhi = 0. ;
	for(int i = 0 ; i < a ; i++)
	{
		double ni = normal() ;
		xhi += ni*ni ;
	}
	return xhi ;
}

double RandomNumber::erlang(double a, int b)
{
	double erlng = 0. ;
	for(int i = 0 ; i < b ; i++)
	{
		double ui = uniform() ;
		erlng += ui*ui ;
	}
	return - a * erlng ;
}


}

































