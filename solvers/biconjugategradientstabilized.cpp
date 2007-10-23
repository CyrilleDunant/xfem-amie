#include "biconjugategradientstabilized.h"
#include "inversediagonal.h"
using namespace Mu ;

BiConjugateGradientStabilized::BiConjugateGradientStabilized(const CoordinateIndexedSparseMatrix &A_, const Vector &b_) :LinearSolver(A_, b_) { }

Vector & BiConjugateGradientStabilized::solve(const Vector &x0, const Preconditionner * precond, const double epsilon , const int maxit , bool verbose )
{
	const Preconditionner * P ;

	
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
	
	Vector r = b - A*x ;
	Vector r_(r) ;
	double rho = std::inner_product(&r[0], &r[r.size()], &r_[0], double(0)) ;
	
	Vector invDiag(r.size()) ;
	
	if(rho < epsilon)
		return x ;
	
	Vector p(r) ;
	Vector p_(r) ;
	P->precondition(p,p_) ;
	//Vector p_ = precondition(p) ;
	
	Vector v = A*p_ ;
	
	double alpha = rho/std::inner_product(&r_[0], &r_[r_.size()], &v[0], double(0)) ;
	
	Vector s = r - v*alpha ;
	
	if(std::abs(s.max()) < epsilon)
	{
		x += p_*alpha ;
		return x ;
	}

	Vector s_(s) ;
	P->precondition(s,s_) ;
	
	Vector t = A*s_ ;
	double omega = std::inner_product(&t[0], &t[t.size()], &s[0], double(0))/std::inner_product(&t[0], &t[t.size()], &t[0], double(0)) ;
	x += p_*alpha +omega*s_ ;
	r = s- t*omega ;
	double rho_ =rho ;
	
	int nit = 0 ;
	int lastit = maxit ;
	if(maxit< 0)
		lastit = b.size()/4 ;
	
	while(nit < lastit)
	{
		nit++ ;
		
		rho = std::inner_product(&r[0], &r[r.size()], &r_[0], double(0)) ;
		if(std::abs(rho) < epsilon)
		{
			if(verbose)
				std::cerr << "\n converged after " << nit << " iterations. Error : " << rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
			
			return x ;
		}
		
		double beta = (rho/rho_)*(alpha/omega) ;
		p = r + (p-v*omega)*beta ;
		
		P->precondition(p,p_) ;
		//p_ = precondition(p) ;
		v = A*p_ ;
		alpha = rho/std::inner_product(&r_[0], &r_[r_.size()], &v[0], double(0)) ;
		s = r - v*alpha ;
		
		//s_ = precondition(s) ;
		P->precondition(s,s_) ;
		t = A*s_ ;
		omega = std::inner_product(&t[0], &t[t.size()], &s[0], double(0))/std::inner_product(&t[0], &t[t.size()], &t[0], double(0)) ;
		
		if(std::abs(omega) < epsilon)
		{
			if(verbose)
				std::cerr << "\n converged after " << nit << " iterations. Error : " << rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
			
			return x ;
		}
		
		x += p_*alpha +s_*omega ;
		r = s- t*omega ;
		rho_ = rho ;

		if(	verbose && nit%100 == 0)
		{
			std::cerr << "\r iteration : " << nit << " error :"<< rho  << "             "<< std::flush ;
		}
		
	}

	if(verbose)
		std::cerr << "\n converged after " << nit << " iterations. Error : " << rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
	
	return x ;
}
