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
*/
class RandomNumber
{
protected:

public:
	RandomNumber() ;

	void reset() { std::srand((unsigned) time(0)) ; } ;
	void reset(int seed) { std::srand(seed) ; } ;

	double uniform() {return (double) std::rand() / (double) RAND_MAX ; } ;
	double uniform(double a) {return a*uniform() ; } ;
	double uniform(double a, double b) {return a + (b-a)*uniform() ; } ;

	int dice(int faces) {return (int) uniform((double) faces) ; } ;

	double triangular() {return uniform() + uniform() ; } ;
	double triangular(double a) {return a * triangular() ; } ;
	double triangular(double a, double b) {return a + (b-a) * triangular() ; } ;

	double exponential(double a) {return -a * std::log(uniform()) ; } ;
	double extreme_value(double a, double b) {return a - b * std::log(- std::log(uniform())) ; } ;

	double normal() ;
	double normal(double a, double b) { return a + sqrt(b) * normal() ; } ;
	double lognormal(double a, double b) { return std::exp(a + b * normal()) ; };

	double logistic(double a, double b) ;

	double cauchy(double a, double b) { return a + b * std::tan(M_PI*(uniform()-0.5)) ; } ;

	double xhi2(int a) ;

	double erlang(double a, int b) ;

	double pareto(double a, double b) { return a * std::pow(1.- uniform(), -1./b) ; } ;

	double weibull(double a, double b) { return a * std::pow(- std::log(uniform()), 1./b) ;};

} ;











}

#endif
