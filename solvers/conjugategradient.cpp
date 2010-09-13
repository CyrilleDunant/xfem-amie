//
// C++ Implementation: conjugategradient
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//



#include "conjugategradient.h"
#include "inversediagonal.h"
#include "tridiagonal.h"
#include "gaussseidellstep.h"
#include "incompletecholeskidecomposition.h"
#include "eigenvalues.h"
#include <limits>

using namespace Mu ;

ConjugateGradient::ConjugateGradient(const CoordinateIndexedSparseMatrix &A_, Vector &b_) :LinearSolver(A_, b_), r(b_.size()),z(b_.size()),p(b_.size()) ,q(b_.size()), cleanup(false), P(NULL), nit(0) { };

bool ConjugateGradient::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
	double realeps = 1e-12 ;
	size_t Maxit ;
	if(maxit != -1)
		Maxit = maxit ;
	else
		Maxit = b.size() ;
	if(x0.size() == b.size())
	{
		x = x0 ;
	}
	else
	{
// 		if(x0.size())
// 		{
// 			std::cout << "ouch" << std::endl ;
// 			exit(0) ;
// 		}
		x = 0 ;
		for(size_t i = 0 ; i < std::min(b.size(), x0.size()) ; i++)
			x[i] = x0[i] ;
	}


	if(precond == NULL && !cleanup)
	{
		cleanup = true ;
// 		P = new InCompleteCholesky(A) ;
// 		P = new InverseDiagonal(A) ;
 		P = new InverseLumpedDiagonal(A) ;
// 		P = new TriDiagonal(A) ;
// 		P = new NullPreconditionner() ;
// 		P = new GaussSeidellStep(A) ;
	}
	else if (precond != NULL)
	{
		delete P ;
		cleanup = false ;
		P = precond ;
	}

	assign(r, A*x-b) ;
	double err0 = std::abs(r).max() ;
	r*=-1 ;

	if (err0 < realeps)
	{
		if(verbose)
			std::cerr << "\n CG "<< p.size() << " converged after " << nit << " iterations. Error : " << err0 << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
		return true ;
	}
	//*************************************
	
	int vsize = r.size() ;
	z = r ;
	P->precondition(r,z) ;

	p = z ;
	q = A*p ;
	
	double last_rho = parallel_inner_product(&r[0], &z[0], vsize) ;
	double alpha = last_rho/parallel_inner_product(&q[0], &p[0], vsize);

	x += p*alpha ;
	r -= q*alpha ;
	nit++ ;
	int n = 0 ;
	//****************************************
	
	assign(r, A*x-b) ;
	r *= -1 ;	
	err0 = std::abs(r).max() ;
	if (err0 < realeps)
	{
		if(verbose)
			std::cerr << "\n CG "<< p.size() << " converged after " << nit << " iterations. Error : " << err0 << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;

		return true ;
	}
	err0 = sqrt( parallel_inner_product(&r[0], &r[0], vsize)) ;
	double neps = std::max(err0*realeps, realeps*realeps) ;
	while(last_rho*last_rho> realeps*realeps && n < Maxit )
	{
		P->precondition(r,z) ;
		
		double rho = parallel_inner_product(&r[0], &z[0], vsize) ;
		double beta = rho/last_rho ;
		p *= beta ;
		p += z ;
		assign(q, A*p) ;
		alpha = rho/parallel_inner_product(&q[0], &p[0], vsize);
		
		r -= q*alpha ;
		x += p*alpha ;
		n++ ;
		if(n%64 == 0)
		{
			assign(r, A*x-b) ;
			r *= -1 ;
			rho =  parallel_inner_product(&r[0], &r[0], vsize) ;
		}
		if(	verbose && nit%64 == 0)
		{
			std::cerr << /*"\r CG "<< p.size() << " iteration : " << nit << " error :"<<*/  sqrt(rho*rho)  /*<< "             "*/<< std::endl /*std::flush*/ ;
		}
// 		if(	verbose )
// 		{
// 			Vector rr = A*x-b ;
// 			for(size_t i = 0 ; i < rr.size() ; i++)
// 			{
// 				std::cerr << rr[i] << "  " << std::flush ;
// 			}
// 			std::cerr << std::endl ;
// 		}
// 		if(	verbose )
// 		{
// 			Vector rr = A*x-b ;
// 			double e = sqrt( parallel_inner_product(&rr[0], &rr[0], vsize)) ;
// 			std::cerr << e << " vs " << sqrt(rho) << std::endl ;
// 		}
		
		last_rho = rho ;
		nit++ ;
		
	}
	assign(r,A*x-b) ;
	double err = sqrt( parallel_inner_product(&r[0], &r[0], vsize)) ;
	
	if(verbose)
	{
		if(nit <= Maxit && last_rho*last_rho< realeps*realeps)
			std::cerr << "\n CG " << p.size() << " converged after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
		else
			std::cerr << "\n CG " << p.size() << " did not converge after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
	}
	
	
	return nit < Maxit && last_rho*last_rho< realeps*realeps;
}

