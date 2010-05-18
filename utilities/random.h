//
// C++ Interface: random generator
//
// Description: simple random generator, using the c++ standard rand() function
//
//
// Author: Alain Giorla
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __RANDOM_H__
#define __RANDOM_H__

#include "../features/features.h"

namespace Mu
{


/** \brief RandomNumber generator
* This class provides a quick access to various random distributions, by 
* using the standard c++ rand() function. When a RandomNumber object is instantiated
* the seed of the pseudo-random suite is reset according to the local time machine.
* In the following description, u refers to a random number between 0 and 1 drawn from a
* uniform distribution, and n to a normal law with parameter 0 and 1.
*/
class RandomNumber
{
protected:

public:
	/**
	* \brief simple constructor
	*/
	RandomNumber() ;

	/**
	* \brief reset the seed based on machine time
	*/
	virtual void reset() { std::srand((unsigned) time(0)) ; } ;
	/**
	* \brief reset the seed according to an arbitrary number
	*/
	virtual void reset(int seed) { std::srand(seed) ; } ;

	/**
	* @return a double between 0 and 1 chosen from an uniform distribution
	*/
	virtual double uniform() {return (double) std::rand() / (double) RAND_MAX ; } ;
	/**
	* @return a double between 0 and a chosen from an uniform distribution
	*/
	virtual double uniform(double a) {return a*uniform() ; } ;
	/**
	* @return a double between a and b chosen from an uniform distribution
	*/
	virtual double uniform(double a, double b) {return a + (b-a)*uniform() ; } ;
	/**
	* @return an integer between 0 and the number of faces (useful for shuffling vectors)
	*/
	virtual int dice(int faces) {return (int) uniform((double) faces) ; } ;

	/**
	* @return a double between 0 and 1, following a triangular probability density function centered on 0.5
	*/
	virtual double triangular() {return uniform() + uniform() ; } ;
	/**
	* @return a double between 0 and a, following a triangular probability density function centered on a/2
	*/
	virtual double triangular(double a) {return a * triangular() ; } ;
	/**
	* @return a double between a and b, following a triangular probability density function centered on (a+b)/2
	*/
	virtual double triangular(double a, double b) {return a + (b-a) * triangular() ; } ;

	/**
	* @param a the rate of the distribution (equals to the inverse of the mean value of the distribution).
	* @return a positive random number following an exponential distribution: -a*ln(u)
	*/
	virtual double exponential(double a) {return -a * std::log(uniform()) ; } ;
	/**
	* @param a the mean of the 
	* @return a random number following an extreme value distribution: a-b*ln(-ln(u))
	*/
	virtual double extreme_value(double a, double b) {return a - b * std::log(- std::log(uniform())) ; } ;

	/**
	* @return a random number following an normal distribution (0,1)
	*/
	virtual double normal() ;
	/**
	* @return a random number following an normal distribution centered on a, width controled by b
	*/
	virtual double normal(double a, double b) { return a + sqrt(b) * normal() ; } ;
	/**
	* @return a random number following a normal distribution (in a logarithmic scale)
	*/
	virtual double lognormal(double a, double b) { return std::exp(a + b * normal()) ; };

	/**
	* @return a random number following an logistic distribution a+b*ln(u/(1-u)), centered on a, width controled by b
	*/
	virtual double logistic(double a, double b) ;
	/**
	* @return a random number following a Cauchy distribution a+b*tan(PI*(u-0.5)) , centered on a, width controled by b
	*/
	virtual double cauchy(double a, double b) { return a + b * std::tan(M_PI*(uniform()-0.5)) ; } ;
	virtual double xhi2(int a) ;
	virtual double erlang(double a, int b) ;
	/**
	* @return a random number following a Pareto distribution a*(1-u^(-1/b))
	*/
	virtual double pareto(double a, double b) { return a * (1.- std::pow(uniform(), -1./b)) ; } ;
	/**
	* @return a random number following a weibull distribution a*(-ln(u)^(1/b))
	*/
	virtual double weibull(double a, double b) { return a * std::pow(- std::log(uniform()), 1./b) ;};

} ;


/**
* Alternative random generator which do NOT rely on the rand() standard function.
*/
class RandomGenerator : public RandomNumber
{
protected:
	long n ;
	long max ;
	long step ;

public:
	RandomGenerator() ;
	RandomGenerator(long seed) ;
	RandomGenerator(long seed, long m, long s) ;

	virtual void reset() ;
	virtual void reset(long seed) ;
	virtual double uniform() ;

	long getMax() {return max ; }
	long getStep() {return step ; }
	void setMax(long m) ;
	void setStep(long s) ;

} ;








}

#endif
