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

using namespace Mu ;

ConjugateGradient::ConjugateGradient(const CoordinateIndexedSparseMatrix &A_, const Vector &b_) :LinearSolver(A_, b_) { };

bool ConjugateGradient::solve(const Vector &x0, const Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
	size_t nit = 0  ;
	size_t Maxit ;
	const Preconditionner * P ;
// 	if(maxit > 0)
// 		Maxit = maxit ;
// 	else
		Maxit = b.size() ;
	
	bool cleanup = false ;
	
	x.resize(b.size(), 0.) ;
	
	if(x0.size() == b.size())
	{
		x = x0 ;
	}
	
	if(precond == NULL)
	{
		cleanup = true ;
// 		P = new InCompleteCholesky(A) ;
		P = new InverseDiagonal(A) ;
// 		P = new NullPreconditionner() ;
	}
	else
		P = precond ;	

	Vector r = b-A*x ;
	double err = std::abs(r).max() ;
	if (err < eps*eps)
	{
		std::cerr << "b in : " << b.min() << ", " << b.max() << std::endl ;
		if(cleanup)
		{
			delete P ;
		}
		return true ;
	}
	//*************************************
	
	Vector z(r) ;
	P->precondition(r,z) ;
	Vector p = z ;
	Vector q = A*p ;
	
	double last_rho = std::inner_product(&r[0], &r[r.size()], &z[0], (double)(0)) ;
	double alpha = last_rho/std::inner_product(&q[0],&q[q.size()] ,&p[0], (double)(0));
	double rho_0 = last_rho ;
	x += p*alpha ;
	r -= q*alpha ;
	err = std::abs(r).max() ;
	
	//****************************************
	
	while(std::abs(last_rho)> rho_0*eps*eps && nit < Maxit)
	{
		P->precondition(r,z) ;
		
		double rho = std::inner_product(&r[0], &r[r.size()], &z[0], 0.) ;
		double beta = rho/last_rho ;
		p = z + p*beta ;
		q = A*p ;
		
		assert(std::inner_product(&q[0],&q[q.size()] ,&p[0], 0.) != 0) ;
		
		alpha = rho/std::inner_product(&q[0],&q[q.size()] ,&p[0], 0.);
		
		r -= q*alpha ;
		x += p*alpha ;
		
		if(	verbose && nit%100 == 0)
		{
			r = b-A*x ;
			std::cerr << "\r iteration : " << nit << " error :"<< last_rho  << "             "<< std::flush ;
		}
		
		last_rho = rho ;
		nit++ ;
		
	}
	r = b-A*x ;
	err = std::abs(r).max() ;
	
	if(verbose)
		std::cerr << "\n converged after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
	
	if(cleanup)
	{
		delete P ;
	}
	
	return (nit < Maxit) ;
}

