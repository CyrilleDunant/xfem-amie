#include "eigenvalues.h"
#include "biconjugategradientstabilized.h"
#include "conjugategradient.h"
#include <numeric>
#include <valarray>

// using namespace Mu ;

namespace Mu
{
double largestEigenValue(const Mu::CoordinateIndexedSparseMatrix & A, bool sym)
  {
		
    std::valarray<double> x_(A.row_size.size()*A.stride) ;
    for(size_t i = 0 ; i< x_.size() ; i++)
    {
      x_[i] = (double)rand()/RAND_MAX ;
    }
    x_/= sqrt(std::inner_product(&x_[0], &x_[x_.size()], &x_[0], double(0)))  ;
    std::valarray<double> x(x_) ;
    
    double eps = 1 ;
    while ( eps > 1e-6 )
    {
      x = A*x_ ;
      x/= sqrt(std::inner_product(&x[0], &x[x.size()], &x[0], double(0)))  ;
      std::valarray<double> delta = x-x_ ;
      eps = sqrt(std::inner_product(&delta[0], &delta[delta.size()], &delta[0], double(0))) ;
      x_ = x ;
    }
    
    x = A*x_ ;
    return std::inner_product(&x[0], &x[x.size()], &x_[0], double(0)) ;
  }

double smallestEigenValue(const Mu::CoordinateIndexedSparseMatrix & A, bool sym)
{
	if(sym)
	{
		Vector x_(A.accumulated_row_size.size()*A.stride) ;
		for(size_t i = 0 ; i< x_.size() ; i++)
		{
			x_[i] = (2.*(double)rand()/RAND_MAX-1) ;
		}
		x_/= sqrt(std::inner_product(&x_[0], &x_[x_.size()], &x_[0], double(0)))  ;
		Vector x(x_) ;
			
		double eps = 1 ;
		size_t it = 0 ;
		while ( eps > 1e-5 && it < 30)
		{
			x = A*x_ ;
			Mu::ConjugateGradient cg(A, x) ;
			cg.solve(x, NULL, 1e-14);
			x = cg.x ;
			x /= sqrt(std::inner_product(&x[0], &x[x.size()], &x[0], double(0)))  ;
			Vector delta = x-x_ ;
			eps = sqrt(std::inner_product(&delta[0], &delta[delta.size()], &delta[0], double(0))) ;
			x_ = x ;
			it++ ;
// 			if(it%100 == 0)
				std::cerr << eps << std::endl ;
		}
			
		x = A*x_ ;
		return std::inner_product(&x[0], &x[x.size()], &x_[0], double(0)) ;
	}
	Vector x_(A.accumulated_row_size.size()*A.stride) ;
	for(size_t i = 0 ; i< x_.size() ; i++)
	{
		x_[i] = (2.*(double)rand()/RAND_MAX-1) ;
	}
	x_/= sqrt(std::inner_product(&x_[0], &x_[x_.size()], &x_[0], double(0)))  ;
	Vector x(x_) ;
    
	double eps = 1 ;
	size_t it = 0 ;
	while ( eps > 1e-5 && it < 30)
	{
		x = A*x_ ;
		Mu::BiConjugateGradientStabilized cg(A, x) ;
		cg.solve(x, NULL, 1e-14);
		x = cg.x ;
		x /= sqrt(std::inner_product(&x[0], &x[x.size()], &x[0], double(0)))  ;
		Vector delta = x-x_ ;
		eps = sqrt(std::inner_product(&delta[0], &delta[delta.size()], &delta[0], double(0))) ;
		x_ = x ;
		it++ ;
		if(it%100 == 0)
			std::cerr << eps << std::endl ;
	}
    
	x = A*x_ ;
	return std::inner_product(&x[0], &x[x.size()], &x_[0], double(0)) ;
}
}