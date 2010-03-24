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

ConjugateGradient::ConjugateGradient(const CoordinateIndexedSparseMatrix &A_, Vector &b_) :LinearSolver(A_, b_) { };

bool ConjugateGradient::solve(const Vector &x0, const Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
	size_t nit = 0  ;
	size_t Maxit ;
	const Preconditionner * P ;
	InverseDiagonal P_alt(A) ;
	if(maxit != -1)
		Maxit = maxit ;
	else
		Maxit = b.size() ;
	
	bool cleanup = false ;
	
	x.resize(b.size(), 0.) ;
	
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


	if(precond == NULL)
	{
		cleanup = true ;
// 		P = new InCompleteCholesky(A) ;
		P = new InverseDiagonal(A) ;
// 		P = new TriDiagonal(A) ;
// 		P = new NullPreconditionner() ;
// 		P = new GaussSeidellStep(A) ;
	}
	else
		P = precond ;	

	Vector r = A*x-b ;
	double err = std::abs(r).max() ;
	r*=-1 ;

	if (err < eps)
	{
		std::cerr << "b in : " << b.min() << ", " << b.max() << ", err = "<< err << std::endl ;
		if(cleanup)
		{
			delete P ;
		}
		return true ;
	}
	//*************************************
	
	int vsize = r.size() ;
	Vector z(r) ;
	P->precondition(r,z) ;
	P_alt.precondition(r,z) ;
	Vector p = z ;
	Vector q = A*p ;
	
	double last_rho = parallel_inner_product(&r[0], &z[0], vsize) ;
	double alpha = last_rho/parallel_inner_product(&q[0], &p[0], vsize);

	x += p*alpha ;
	r -= q*alpha ;
	
	//****************************************
	double neps = 1e-9 ;
	while(std::abs(last_rho)> std::max(err*neps*neps, neps*neps) && nit < Maxit )
	{
		P->precondition(r,z) ;
		P_alt.precondition(r,z) ;
		
		double rho = parallel_inner_product(&r[0], &z[0], vsize) ;
		double beta = rho/last_rho ;
		p *=beta ;
		p += z ;
		assign(q, A*p) ;
		alpha = rho/parallel_inner_product(&q[0], &p[0], vsize);
		
		r -= q*alpha ;
		x += p*alpha ;
		
		if(nit%64 == 0)
		{
			assign(r, A*x-b) ;
			r *= -1 ;		
		}
		if(	verbose && nit%64 == 0)
		{
			std::cerr << "\r iteration : " << nit << " error :"<<  std::abs(rho)  << "             "<< std::flush ;
		}
		
		last_rho = rho ;
		nit++ ;
		
	}
	assign(r,A*x-b) ;
	err = std::abs(r).max() ;
	
	if(verbose)
	{
		if(nit < Maxit)
			std::cerr << "\n converged after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
		else
			std::cerr << "\n did not converge after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
	}
	
	if(cleanup)
	{
		delete P ;
	}
	
	return (err < eps)  ;
}

