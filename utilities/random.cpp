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






RandomGenerator::RandomGenerator()
{
	n = 1 ;
	max = (long) (std::pow(2,31) - 1) ;
	step = (long) (std::pow(7,5)) ;
	reset() ;
}

RandomGenerator::RandomGenerator(long seed)
{
	n = seed ;
	max = (long) (std::pow(2,31) - 1) ;
	step = (long) (std::pow(7,5)) ;
}

RandomGenerator::RandomGenerator(long seed, long m, long s)
{
	n = seed ;
	setMax(m) ;
	setStep(s) ;
}

void RandomGenerator::reset()
{
	n = time(0) ;
	if(n < 1)
		n = 1 ; 
}

void RandomGenerator::reset(long seed)
{
	if(seed > 0)
		n = seed ;
}

void RandomGenerator::setMax(long m)
{
	if(m > 0)
		max = m ;
}

void RandomGenerator::setStep(long s)
{
	if(s > 0 && s%max != 0)
		step = s ;
}

double RandomGenerator::uniform()
{
	n *= step ;
	n = n%max ;
	if(n==0)
		n = 1 ;
	return (double) (n/max) ;
}


RandomDistribution::RandomDistribution()
{
	rnd = new RandomNumber() ;
}

UniformDistribution::UniformDistribution() : RandomDistribution()
{
	a = 0. ;
	b = 1. ;
}
UniformDistribution::UniformDistribution(double x) : RandomDistribution()
{
	a = 0. ;
	b = x ;
}
UniformDistribution::UniformDistribution(double x, double y) : RandomDistribution()
{
	a = x ;
	b = y ;
}

TriangularDistribution::TriangularDistribution() : UniformDistribution() 
{
}
TriangularDistribution::TriangularDistribution(double x) : UniformDistribution(x)
{
}
TriangularDistribution::TriangularDistribution(double x, double y) : UniformDistribution(x,y)
{
}

ExponentialDistribution::ExponentialDistribution(double x) : RandomDistribution()
{
      a = x ;
}

ExtremeValueDistribution::ExtremeValueDistribution(double x, double y) : NormalDistribution(x,y)
{
}

NormalDistribution::NormalDistribution() : UniformDistribution()
{
}
NormalDistribution::NormalDistribution(double x, double y) : UniformDistribution(x,y)
{
}

LogNormalDistribution::LogNormalDistribution() : NormalDistribution()
{
}
LogNormalDistribution::LogNormalDistribution(double x, double y) : NormalDistribution(x,y)
{
}

LogisticDistribution::LogisticDistribution(double x, double y) : NormalDistribution(x,y)
{
}

CauchyDistribution::CauchyDistribution(double x, double y) : NormalDistribution(x,y)
{
}

Xhi2Distribution::Xhi2Distribution(int x) : RandomDistribution()
{
      n = x ;
}

ErlangDistribution::ErlangDistribution(double x, int y) : Xhi2Distribution(y)
{
      a = x ;
}

ParetoDistribution::ParetoDistribution(double x, double y) : NormalDistribution(x,y)
{
}

WeibullDistribution::WeibullDistribution(double x, double y) : NormalDistribution(x,y)
{
}



}






























