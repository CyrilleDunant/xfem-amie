//
// C++ Implementation: polakribiereconjugategradient
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "polakribiereconjugategradient.h"
#include "inversediagonal.h"
#include <limits>

using namespace Amie ;

ConjugateGradientWithSecant::ConjugateGradientWithSecant(Assembly * a, size_t n) : NonLinearSolver(a), nssor(n) { } 


bool ConjugateGradientWithSecant::solve(const Vector &x0, Preconditionner * precond , const double eps , const int maxit , bool verbose )
{
	double sigma = 4;//(1.+sqrt(5.))/2. ; //secant parameter
	size_t resetParameter = 16 ;
	size_t maxSecantIteration = 16 ;
	
	size_t k = 0 ;
	size_t i = 0 ;
	const CoordinateIndexedSparseMatrix &A = assembly->getMatrix() ;
	Vector & b = assembly->getForces() ;
	Vector & y = assembly->getDisplacements();
	bool nl = assembly->nonLinearStep() ;

	
	if(!nl) 
	{
		std::cerr << "Linear problem, falling back to linear CG" << std::endl;
		ConjugateGradient cg(assembly) ;
		cg.colstart = this->colstart ;
		cg.rowstart = this->rowstart ;
                cg.nssor = this->nssor ;
		bool ret =  cg.solve(x0, nullptr,eps, -1, verbose) ;
		x.resize(cg.x.size())  ;
		x = cg.x ;
		return ret ;
	}
	
	x.resize(b.size()) ;
	if(x0.size() == b.size())
	{
		x = x0 ;
	}
	else
	{
		x = 0 ;
		for(size_t i = 0 ; i < std::min(b.size(), x0.size()) ; i++)
			x[i] = x0[i] ;
	}

	InverseDiagonal P(A) ;
// 	InverseLumpedDiagonal P(A) ;
// 	InCompleteCholesky P(A) ;
// 	NullPreconditionner P ;
	
	size_t vsize = b.size() ;
	Vector & nlb = assembly->getNonLinearForces() ;
	Vector r = ((A+ assembly->getNonLinearMatrix())*y)-b-nlb;
	Vector s(r);
	P.precondition(r,s) ;//ConjugateGradient(A,r).solve(x, nullptr, 1e-12, 200, true) ;// 
// 	s*=-1 ;
	Vector d = s ;
	
	double delta_new = parallel_inner_product(&r[0],&d[0], vsize );
	double delta_0 = delta_new ;
	
	std::cerr << " iteration : " << i << " error :"<< delta_new <<", max : "  << x.max() << ", min : "  << x.min() << "             "<< std::endl ;
	
	
	while((i < b.size()/4) && 
	      (delta_new > eps*eps*delta_0) && 
	      (delta_new > 4.*POINT_TOLERANCE*POINT_TOLERANCE))
	{
		bool successefulSecant = false ;
		size_t secantcount = 0;
		
		while(!successefulSecant && secantcount<1)
		{
			assembly->nonLinearStep() ;
			size_t j = 0 ;
			double delta_d = parallel_inner_product(&d[0],&d[0], vsize );
			double alpha = -sigma ;
			Vector temp = -(Vector)((A+ assembly->getNonLinearMatrix())*(y + d*sigma)-nlb-b) ;
			
			double eta_prev = parallel_inner_product(&temp[0], &d[0], vsize) ;
			
			do
			{
				temp = -(Vector)((A+ assembly->getNonLinearMatrix())*y-nlb-b) ;
				double eta =  parallel_inner_product(&temp[0], &d[0], vsize) ;
				
				if(std::abs(eta_prev-eta) < std::numeric_limits<double>::epsilon())
				{
					sigma *=4 ;
					secantcount++ ;
					break ;
				}
				
				
				alpha *= eta/(eta_prev-eta) ;
				y += d*alpha ;
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
			
			std::cerr << "eta_prev = " << eta_prev << std::endl ;
			secantcount++ ;
		}
		
		r = (A+assembly->getNonLinearMatrix())*y -b-nlb;
		
		double delta_old = delta_new ;
		double delta_mid = parallel_inner_product(&r[0], &s[0], vsize) ;
		
		P.precondition(r,s) ;//s = ConjugateGradient(A+assembly->getNonLinearMatrix(), r).solve(s, &P, 1e-24, 10, false) ;//
// 		s*=-1 ;
		delta_new = parallel_inner_product(&r[0], &s[0], vsize) ;
		std::cerr << i << " :: delta_new = " << delta_new << std::endl ;
		std::cerr << i << " :: delta_mid = " << delta_mid << std::endl ;
		std::cerr << i << " :: delta_old = " << delta_old << std::endl ;
		
		double beta = (delta_new-delta_mid)/delta_old ;
		k++ ;
		
		if(i%128 == 0)
		{
			
			std::cerr << /*"\r iteration : " <<*/ i << "  "<< delta_new  << "                       "<< std::endl ;
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
	
	r = (A+assembly->getNonLinearMatrix())*y -b-nlb;
	
	
	std::cerr << "\n iteration : " << i << " error :"<< std::abs(r).max() <<", max : "  << y.max() << ", min : "  << y.min() << "             "<< std::endl ;
	
	x.resize(assembly->getDisplacements().size()) ;
	x = assembly->getDisplacements() ;
	return i < b.size()/4;
}

