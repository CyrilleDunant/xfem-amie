#include "biconjugategradientstabilized.h"
#include "inversediagonal.h"
using namespace Mu ;

BiConjugateGradientStabilized::BiConjugateGradientStabilized(const CoordinateIndexedSparseMatrix &A_, Vector &b_) :LinearSolver(A_, b_) { }

bool BiConjugateGradientStabilized::solve(const Vector &x0, const Preconditionner * precond, const double epsilon , const int maxit , bool verbose )
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
	
	Vector r = A*x-b ;
	r *= -1 ;
	Vector r_(r) ;
	P->precondition(r,r_) ;
	double rho = std::inner_product(&r[0], &r[r.size()], &r_[0], double(0)) ;
	
	Vector invDiag(r.size()) ;
	
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
	double omega = parallel_inner_product(&t[0], &s[0], vsize)/parallel_inner_product(&t[0], &t[0], vsize) ;
	x += p_*alpha +omega*s_ ;
	r = s- t*omega ;
	double rho_ =rho ;
	
	int nit = 0 ;
	int lastit = std::min(maxit, (int)b.size()/4) ;
	if(maxit< 0)
		lastit = b.size() ;
	
	while(nit < lastit)
	{
		nit++ ;
		
		rho = parallel_inner_product(&r[0], &r_[0], vsize) ;
		if(std::abs(rho) < epsilon*epsilon)
		{
			if(verbose)
				std::cerr << "\n converged after " << nit << " iterations. Error : " << rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
			if(cleanup)
				delete P ;
			
			return true ;
		}
		
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
		omega = parallel_inner_product(&t[0], &s[0], vsize)/parallel_inner_product(&t[0], &t[0], vsize) ;
		
		if(std::abs(omega) < epsilon*epsilon)
		{
			if(verbose)
				std::cerr << "\n converged after " << nit << " iterations. Error : " << rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
			if(cleanup)
				delete P ;
			return true ;
		}
		
		x += p_*alpha +s_*omega ;
		r = s- t*omega ;
		rho_ = rho ;

		if(verbose && nit%60 == 0)
		{
// 			r = b - A*x ;
			std::cerr << "\r iteration : " << nit << " error :"<< rho  << "             "<< std::flush ;
		}
		
	}

	if(verbose)
		std::cerr << "\n did not converge after " << nit << " iterations. Error : " << rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
	
	if(cleanup)
		delete P ;

	return false ;
}
