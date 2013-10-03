//
// C++ Implementation: conjugategradient
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
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
#include <sys/time.h>

using namespace Mu ;

ConjugateGradient::ConjugateGradient(const CoordinateIndexedSparseMatrix &A_, Vector &b_) :LinearSolver(A_, b_), r(b_.size()),z(b_.size()),p(b_.size()) ,q(b_.size()), cleanup(false), P(nullptr), nit(0) { };

bool ConjugateGradient::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
	double realeps = 1e-12 ;
	size_t Maxit ;
	if(maxit != -1)
		Maxit = maxit ;
	else
		Maxit = std::max(round(b.size()), 1000.) ;
	if(x0.size() == b.size())
	{
		x.resize(x0.size());
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

	if(precond == nullptr && !cleanup)
	{
		cleanup = true ;
// 		P = new InCompleteCholesky(A) ;
		P = new InverseDiagonal(A) ;
//  		P = new InverseLumpedDiagonal(A) ;
// 		P = new TriDiagonal(A) ;
// 		P = new NullPreconditionner() ;
// 		P = new GaussSeidellStep(A) ;
	}
	else if (precond != nullptr)
	{
		delete P ;
		cleanup = false ;
		P = precond ;
	}

	assign(r, A*x-b, rowstart, colstart) ;
	int vsize = r.size() ;
	double err0 = sqrt( parallel_inner_product(&r[0], &r[0], vsize)) ;
	r*=-1 ;

	if (err0 < eps)
	{
		if(verbose)
			std::cerr << "\n CG "<< p.size() << " converged after " << nit << " iterations. Error : " << err0 << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
		return true ;
	}
	//*************************************
	
	
	z = r ;
	P->precondition(r,z) ;

	p = z ;
	q = A*p ;
	
	double last_rho = parallel_inner_product(&r[0], &z[0], vsize) ;
	double pq = parallel_inner_product(&q[0], &p[0], vsize);
	double alpha = last_rho/pq ;
	
	x += p*alpha ;
	r -= q*alpha ;
	nit++ ;
	//****************************************
	
	assign(r, A*x-b, rowstart, colstart) ;
	r *= -1 ;
	err0 = sqrt( parallel_inner_product(&r[0], &r[0], vsize)) ;
	if (err0 < eps)
	{
		if(verbose)
			std::cerr << "\n CG "<< p.size() << " converged after " << nit << " iterations. Error : " << err0 << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;

		return true ;
	}
	timeval time0, time1 ;
	gettimeofday(&time0, nullptr);
// 	double neps = /*std::min(*/realeps*realeps/*, err0*realeps)*/ ; //std::max(err0*realeps, realeps*realeps) ;
	double neps = std::max(realeps, eps) ;
	double rho = 0 ;
	double beta = 0 ;
	while(last_rho*last_rho*vsize*vsize > std::max(neps*neps*err0, neps*neps) && nit < Maxit )
	{
		P->precondition(r,z) ;

		rho = parallel_inner_product(&r[0], &z[0], vsize) ;

		beta = rho/last_rho ;
		
		#pragma omp parallel for schedule(runtime) if (vsize > 10000)
		for(size_t i = 0 ; i < vsize ; i++)
			p[i] = p[i]*beta+z[i] ;

		assign(q, A*p, rowstart, colstart) ;
		pq =  parallel_inner_product(&q[0], &p[0], vsize);
		alpha = rho/pq;
		if(std::abs(pq) < POINT_TOLERANCE_2D*POINT_TOLERANCE_2D)
		{
			last_rho = 0 ;
			break ;
		}
		r -= q*alpha ;
		x += p*alpha ;
		
		if(nit%32 == 0)
		{
			assign(r, A*x-b, rowstart, colstart) ;
			r *= -1 ;
		}
		if(	verbose && nit%128 == 0)
		{
			std::cerr << sqrt(rho) << std::endl  ;
		}

		last_rho = rho ;
		nit++ ;
	}
	gettimeofday(&time1, nullptr);
	double delta = time1.tv_sec*1000000 - time0.tv_sec*1000000 + time1.tv_usec - time0.tv_usec ;
	std::cerr << "mflops: "<< nit*((2.+2./32.)*A.array.size()+(4+1./32.)*p.size())/delta << std::endl ;

	assign(r,A*x-b, rowstart, colstart) ;
	double err = sqrt( parallel_inner_product(&r[0], &r[0], vsize)) ;
	
	if(verbose)
	{
		if(nit <= Maxit && last_rho*last_rho< std::max(neps*neps*err0, neps*neps))
			std::cerr << "\n CG " << p.size() << " converged after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
		else
			std::cerr << "\n CG " << p.size() << " did not converge after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
	}
	for (size_t i = 0 ; i < x.size() ; i++)
	{
		std::cout << x[i] << "  " << std::flush ;
	}
	return nit <= Maxit && last_rho*last_rho< std::max(neps*neps*err0, neps*neps);
}

