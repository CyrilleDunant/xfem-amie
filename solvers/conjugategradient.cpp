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
#include <omp.h>

using namespace Amie ;

ConjugateGradient::ConjugateGradient(const CoordinateIndexedSparseMatrix &A_, Vector &b_) :LinearSolver(A_, b_), r(b_.size()),z(b_.size()),p(b_.size()) ,q(b_.size()), cleanup(false), P(nullptr), nit(0) { };

bool ConjugateGradient::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
	double realeps = std::max(1e-9, eps) ;
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
	
	if(rowstart)
	{
		for(size_t i = 0 ; i < rowstart ; i++)
			x[i] = b[i] ;
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
	double err0 = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
	r*=-1 ;


	if (err0 < realeps)
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
	
	double last_rho = parallel_inner_product_restricted(&r[rowstart], &z[rowstart], vsize-rowstart) ;
	double pq = parallel_inner_product_restricted(&q[rowstart], &p[rowstart], vsize-rowstart);
	double alpha = last_rho/pq ;
	
	#pragma omp parallel for schedule(static) if (vsize > 10000)
	for(size_t i = rowstart ; i < vsize ; i++)
	{
		r[i] -= q[i]*alpha ;
		x[i] += p[i]*alpha ;
	}
	nit++ ;
	//****************************************
	
	assign(r, A*x-b, rowstart, colstart) ;
	#pragma omp parallel for schedule(static) if (vsize > 10000)
	for(size_t i = rowstart ; i < vsize ; i++)
			r[i] *= -1 ;

	err0 = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
	if (err0 < realeps)
	{
		if(verbose)
			std::cerr << "\n CG "<< p.size() << " converged after " << nit << " iterations. Error : " << err0 << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;

		return true ;
	}
#ifdef HAVE_OMP
	double t0 = omp_get_wtime() ;
#else
	double t0 = 0 ;
#endif
// 	double neps = /*std::min(*/realeps*realeps/*, err0*realeps)*/ ; //std::max(err0*realeps, realeps*realeps) ;
	double rho = 0 ;
	double beta = 0 ;
	while(last_rho*last_rho*vsize*vsize > std::max(realeps*realeps*err0, realeps*realeps) && nit < Maxit )
	{
		P->precondition(r,z) ;

		rho = parallel_inner_product_restricted(&r[rowstart], &z[rowstart], vsize-rowstart) ;

		beta = rho/last_rho ;

		
		#pragma omp parallel for schedule(static) if (vsize > 10000)
		for(size_t i = rowstart ; i < vsize ; i++)
			p[i] = p[i]*beta+z[i] ;

		assign(q, A*p, rowstart, colstart) ;
		pq =  parallel_inner_product_restricted(&q[rowstart], &p[rowstart], vsize-rowstart);
		alpha = rho/pq;
	
		if(std::abs(pq) < POINT_TOLERANCE_2D*POINT_TOLERANCE_2D)
		{
			last_rho = 0 ;
			break ;
		}
		#pragma omp parallel for schedule(static) if (vsize > 10000)
		for(size_t i = rowstart ; i < vsize ; i++)
		{
			r[i] -= q[i]*alpha ;
			x[i] += p[i]*alpha ;
		}
		
		if(nit%256 == 0)
		{
			assign(r, A*x-b, rowstart, rowstart) ;
			#pragma omp parallel for schedule(static) if (vsize > 10000)
			for(size_t i = rowstart ; i < vsize ; i++)
				r[i] *= -1 ;
		}
		if( verbose && nit%256 == 0 )
		{
			std::cerr << sqrt(rho) << std::endl  ;
		}


		last_rho = rho ;
		nit++ ;
	}
#ifdef HAVE_OMP
	double delta = std::max((omp_get_wtime() - t0)*1e6, 1e-14) ;
#else
	double delta = 1 ;
#endif
	
	std::cerr << "mflops: "<< nit*((2.+2./256.)*A.array.size()+(4+1./256.)*p.size())/delta << std::endl ;

	assign(r,A*x-b, rowstart, rowstart) ;
	double err = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
	
	if(verbose)
	{
		if(nit <= Maxit && last_rho*last_rho< std::max(realeps*realeps*err0, realeps*realeps))
			std::cerr << "\n CG " << p.size() << " converged after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
		else
			std::cerr << "\n CG " << p.size() << " did not converge after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
	}

	return nit <= Maxit && last_rho*last_rho< std::max(realeps*realeps*err0, realeps*realeps);
}

