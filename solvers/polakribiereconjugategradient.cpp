//
// C++ Implementation: polakribiereconjugategradient
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "polakribiereconjugategradient.h"
#include "inversediagonal.h"

using namespace Mu ;

ConjugateGradientWithSecant::ConjugateGradientWithSecant(Assembly * a) : NonLinearSolver(a) { } 


Vector & ConjugateGradientWithSecant::solve(const Vector &x0, const Preconditionner * precond , const double eps , const int maxit , bool verbose )
{
	double sigma = 4;//(1.+sqrt(5.))/2. ; //secant parameter
	size_t resetParameter = 16 ;
	size_t maxSecantIteration = 16 ;
	
	size_t k = 0 ;
	size_t i = 0 ;
	const CoordinateIndexedSparseMatrix &A = assembly->getMatrix() ;
	const Vector & b = assembly->getForces() ;
	Vector & x = assembly->getDisplacements();
	bool nl = assembly->nonLinearStep() ;
	
	if(x.size() != b.size())
	{
		x.resize(b.size()) ;
		x = 0 ;
	}
	
	if(!nl) 
	{
		std::cout << "Linear problem, falling back to linear CG" << std::endl;
		x =  ConjugateGradient(assembly->getMatrix(), b).solve(x, NULL,eps, b.size()/8, true) ;
		return x ;
		
	}
// 	else
// 	{
// 		x =  ConjugateGradient(assembly->getMatrix(), b).solve(x, NULL,eps, b.size(), true) ;
// 		for(size_t i = 0 ; i < 2 ; i++)
// 		{
// 			assembly->nonLinearStep() ;
// 			x = ConjugateGradient(assembly->getMatrix()+assembly->getNonLinearMatrix(), b+assembly->getNonLinearForces()).solve(x, NULL,eps, b.size(), true) ;
// 		}
// 		return x ;
// 		
// 	}
	InverseDiagonal P(A) ;
// 	InverseLumpedDiagonal P(A) ;
// 	InCompleteCholesky P(A) ;
// 	NullPreconditionner P ;
	
	size_t vsize = b.size() ;
	Vector & nlb = assembly->getNonLinearForces() ;
	Vector r = ((A+ assembly->getNonLinearMatrix())*x)-b-nlb;
	Vector s(r);
	P.precondition(r,s) ;//ConjugateGradient(A,r).solve(x, NULL, 1e-12, 200, true) ;// 
// 	s*=-1 ;
	Vector d = s ;
	
	double delta_new = std::inner_product(&r[0],&r[vsize],&d[0], double(0) );
	double delta_0 = delta_new ;
	
	std::cout << " iteration : " << i << " error :"<< delta_new <<", max : "  << x.max() << ", min : "  << x.min() << "             "<< std::endl ;
	
	
	while((i < 800 /*b.size()/4*/) && 
	      (delta_new > eps*eps*delta_0) && 
	      (delta_new > 4.*POINT_TOLERANCE*POINT_TOLERANCE))
	{
		bool successefulSecant = false ;
		size_t secantcount = 0;
		
		while(!successefulSecant && secantcount<1)
		{
			assembly->nonLinearStep() ;
			size_t j = 0 ;
			double delta_d = std::inner_product(&d[0],&d[vsize],&d[0], double(0) );
			double alpha = -sigma ;
			Vector temp = -(Vector)((A+ assembly->getNonLinearMatrix())*(x + d*sigma)-nlb-b) ;
			
			double eta_prev = std::inner_product(&temp[0], &temp[vsize], &d[0], double(0)) ;
			
			do
			{
				Vector temp = -(Vector)((A+ assembly->getNonLinearMatrix())*x-nlb-b) ;
				double eta =  std::inner_product(&temp[0], &temp[vsize], &d[0], double(0)) ;
				
				if(std::isnan(eta/(eta_prev-eta)))
				{
					sigma *=4 ;
					secantcount++ ;
					break ;
				}
				
				
				alpha *= eta/(eta_prev-eta) ;
				x += d*alpha ;
				eta_prev = eta ;
				j++ ;
				
				
				if( (alpha*alpha*delta_d < eps*eps) )
				{
					sigma = 4 ;
					successefulSecant = true ;
					break ;
				}
				
				if(j >= maxSecantIteration)
				{
					sigma *=4 ;
					secantcount++ ;
					break ;
				}
				assembly->nonLinearStep() ;
				
			} while (true) ;
			
			std::cout << "eta_prev = " << eta_prev << std::endl ;
			secantcount++ ;
		}
		
		r = (A+assembly->getNonLinearMatrix())*x -b-nlb;
		
		double delta_old = delta_new ;
		double delta_mid = std::inner_product(&r[0], &r[vsize], &s[0], double(0)) ;
		
		P.precondition(r,s) ;//s = ConjugateGradient(A+assembly->getNonLinearMatrix(), r).solve(s, &P, 1e-24, 10, false) ;//
// 		s*=-1 ;
		delta_new = std::inner_product(&r[0], &r[vsize], &s[0], double(0)) ;
		std::cout << i << " :: delta_new = " << delta_new << std::endl ;
		std::cout << i << " :: delta_mid = " << delta_mid << std::endl ;
		std::cout << i << " :: delta_old = " << delta_old << std::endl ;
		
		double beta = (delta_new-delta_mid)/delta_old ;
		k++ ;
		
		if(i%128 == 0)
		{
			
			std::cout << /*"\r iteration : " <<*/ i << "  "<< delta_new  << "                       "<< std::endl ;
		}
		if(k == resetParameter || beta < 0)
		{
			d = s ;
			k = 0 ;
		}
		else
			d = s+d*beta ;
		
		i++ ;
		
	}
	
	r = (A+assembly->getNonLinearMatrix())*x -b-nlb;
	
	std::cout << "\n iteration : " << i << " error :"<< std::abs(r).max() <<", max : "  << x.max() << ", min : "  << x.min() << "             "<< std::endl ;
	
	return assembly->getDisplacements() ;
}

