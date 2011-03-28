// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011

#include "biconjugategradientstabilized.h"
#include "inversediagonal.h"

#include <sys/time.h>

using namespace Mu ;

BiConjugateGradientStabilized::BiConjugateGradientStabilized(const CoordinateIndexedSparseMatrix &A_, Vector &b_) :LinearSolver(A_, b_) { }

bool BiConjugateGradientStabilized::solve(const Vector &x0, Preconditionner * precond, const double epsilon , const int maxit , bool verbose )
{
	Preconditionner * P ;

	
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
	
	Vector r = A*x-b ;
	r *= -1 ;
	Vector r_(r) ;
	P->precondition(r,r_) ;
	double rho = std::inner_product(&r[0], &r[r.size()], &r_[0], double(0)) ;
	
	
	if(std::abs(rho) < epsilon*epsilon)
		return true ;
	
	Vector p(r) ;
	Vector p_(r) ;
	P->precondition(p,p_) ;
	//Vector p_ = precondition(p) ;
	int vsize = r.size() ;
	Vector v = A*p_ ;
	
	double alpha = rho/parallel_inner_product(&r_[0], &v[0], vsize) ;
	
	Vector s = r - v*alpha ;
	
	if(std::abs(s.max()) < epsilon)
	{
		x += p_*alpha ;
		if(cleanup)
			delete P ;
		return true ;
	}

	Vector s_(s) ;
	
	P->precondition(s,s_) ;
	
	Vector t = A*s_ ;
	Vector t__(t) ;
	Vector s__(s) ;
	P->precondition(t,t__) ;
	P->precondition(s,s__) ;
	double omega = parallel_inner_product(&t__[0], &s__[0], vsize)/parallel_inner_product(&t__[0], &t__[0], vsize) ;
	x += p_*alpha +omega*s_ ;
	r = s- t*omega ;
	double rho_ =rho ;
	
	double err0 = sqrt( std::abs(parallel_inner_product(&r[0], &r[0], vsize))) ;
	
	int nit = 0 ;
	int lastit = std::min(maxit, (int)b.size()/4) ;
	if(maxit< 0)
		lastit = b.size() ;
	timeval time0, time1 ;
	gettimeofday(&time0, NULL);
	
	while(nit < lastit && std::abs(rho)*vsize*vsize > std::max(std::abs(err0)*epsilon*epsilon, epsilon*epsilon) )
	{
		nit++ ;
		
		rho = parallel_inner_product(&r[0], &r_[0], vsize) ;
		
		double beta = (rho/rho_)*(alpha/omega) ;
		p = r + (p-v*omega)*beta ;
		
		P->precondition(p,p_) ;
		//p_ = precondition(p) ;
		assign(v, A*p_) ;
		alpha = rho/parallel_inner_product(&r_[0], &v[0], vsize) ;
		s = r - v*alpha ;
		
		//s_ = precondition(s) ;
		P->precondition(s,s_) ;
		assign(t, A*s_) ;

		if(false) //two variants for preconditionning BiCGSTAB
		{
			omega = parallel_inner_product(&t[0], &s[0], vsize)/parallel_inner_product(&t[0], &t[0], vsize) ;
		}
		else
		{
			P->precondition(t,t__) ;
			P->precondition(s,s__) ;
			omega = parallel_inner_product(&t__[0], &s__[0], vsize)/parallel_inner_product(&t__[0], &t__[0], vsize) ;
		}
		
		x += p_*alpha +s_*omega ;
		r = s- t*omega ;
		rho_ = rho ;

		
		if(verbose && nit%128 == 0)
		{
			std::cerr <<  sqrt(std::abs(rho)) << std::endl  ;
		}
		
	}
	gettimeofday(&time1, NULL);
	double delta = time1.tv_sec*1000000 - time0.tv_sec*1000000 + time1.tv_usec - time0.tv_usec ;
	std::cerr << "mflops: "<< nit*2*((2.)*A.array.size()+6*p.size())/delta << std::endl ;

	assign(r,A*x-b) ;
	double err = sqrt( parallel_inner_product(&r[0], &r[0], vsize)) ;
	
	if(verbose)
	{
		if(nit <= lastit && std::abs(rho) <= std::max(std::abs(err0)*epsilon*epsilon, epsilon*epsilon))
			std::cerr << "\n BiCGStab " << p.size() << " converged after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
		else
			std::cerr << "\n BiCGStab " << p.size() << " did not converge after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
	}
	
	if(cleanup)
		delete P ;
	
	return nit < lastit && std::abs(rho) <= std::max(std::abs(err0)*epsilon*epsilon, epsilon*epsilon);
}
