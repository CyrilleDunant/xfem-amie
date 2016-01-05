// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013

#include "eigenvalues.h"
#include "biconjugategradientstabilized.h"
#include "conjugategradient.h"
#include <numeric>
#include <valarray>


namespace Amie
{
double largestEigenValue(Assembly * a, bool sym)
{
	srand(0) ;
	
	std::valarray<double> x_(a->getMatrix().row_size.size()*a->getMatrix().stride) ;
	for(size_t i = 0 ; i< x_.size() ; i++)
	{
		x_[i] = (double)rand()/RAND_MAX ;
	}
	x_/= sqrt(std::inner_product(&x_[0], &x_[x_.size()], &x_[0], double(0)))  ;
	std::valarray<double> x(x_) ;
	
	double eps = 1. ;
	while ( eps > 1e-6 )
	{
		x = a->getMatrix()*x_ ;
		x/= sqrt(std::inner_product(&x[0], &x[x.size()], &x[0], double(0)))  ;
		std::valarray<double> delta = x-x_ ;
		eps = sqrt(std::inner_product(&delta[0], &delta[delta.size()], &delta[0], double(0))) ;
		x_ = x ;
	}
	
	x = a->getMatrix()*x_ ;
	return std::inner_product(&x[0], &x[x.size()], &x_[0], double(0)) ;
}

double smallestEigenValue(Assembly * a, bool sym)
{
	srand(0) ;
	if(sym)
	{
		Vector x_(a->getMatrix().accumulated_row_size.size()*a->getMatrix().stride) ;
		for(size_t i = 0 ; i< x_.size() ; i++)
		{
			x_[i] = (2.*(double)rand()/RAND_MAX-1) ;
		}
		x_/= sqrt(std::inner_product(&x_[0], &x_[x_.size()], &x_[0], double(0)))  ;
		Vector x(x_) ;
			
		double eps = 1. ;
		size_t it = 0 ;
                Assembly a_(*a) ;
		while ( eps > 1e-6 )
		{
			a_.getForces() = a->getMatrix()*x_ ;
			Amie::ConjugateGradient cg(a) ;
			cg.solve(x, nullptr, 1e-14);
			x = cg.x ;
			x /= sqrt(std::inner_product(&x[0], &x[x.size()], &x[0], double(0)))  ;
			Vector delta = x-x_ ;
			eps = sqrt(std::inner_product(&delta[0], &delta[delta.size()], &delta[0], double(0))) ;
			x_ = x ;
			it++ ;
// 			if(it%100 == 0)
				std::cerr << eps << std::endl ;
		}
			
		x = a->getMatrix()*x_ ;
		return std::inner_product(&x[0], &x[x.size()], &x_[0], double(0)) ;
	}
	Vector x_(a->getMatrix().accumulated_row_size.size()*a->getMatrix().stride) ;
	for(size_t i = 0 ; i< x_.size() ; i++)
	{
		x_[i] = (2.*(double)rand()/RAND_MAX-1) ;
	}
	x_/= sqrt(std::inner_product(&x_[0], &x_[x_.size()], &x_[0], double(0)))  ;
	Vector x(x_) ;
    
	double eps = 1 ;
	size_t it = 0 ;
        Assembly a_(*a) ;
	while ( eps > 1e-5 && it < 30)
	{
		x = a->getMatrix()*x_ ;
                a_.getForces() = x ;
		Amie::BiConjugateGradientStabilized cg(&a_) ;
		cg.solve(x, nullptr, 1e-14);
		x = cg.x ;
		x /= sqrt(std::inner_product(&x[0], &x[x.size()], &x[0], double(0)))  ;
		Vector delta = x-x_ ;
		eps = sqrt(std::inner_product(&delta[0], &delta[delta.size()], &delta[0], double(0))) ;
		x_ = x ;
		it++ ;
		if(it%100 == 0)
			std::cerr << eps << std::endl ;
	}
    
	x = a->getMatrix()*x_ ;
	return std::inner_product(&x[0], &x[x.size()], &x_[0], double(0)) ;
}
}
