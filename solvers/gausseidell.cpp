//
// C++ Implementation: gausseidell
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "gausseidell.h"

namespace Mu {

GaussSeidel::GaussSeidel(const CoordinateIndexedSparseMatrix &A_, const Vector &b_) :LinearSolver(A_, b_) { };

Vector & GaussSeidel::solve(const Vector &x0, const Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
	double err=10 ;
	
	if(x0.size() == b.size())
	{
		x = x0 ;
	}
	
	size_t nit=0 ;
	size_t Maxit ;
	if(maxit > 0)
		Maxit = maxit ;
	else
		Maxit = 2*b.size() ;
	
// 	double omega = 1 ;
	
	Vector inverseDiagonal = A.inverseDiagonal() ;
	
	while((err > eps) && nit<Maxit)
	{
		err = 0 ;
		for(size_t i = 0 ; i < x.size() ; i++)
		{
			double delta = 0 ;
			size_t start_index = A.accumulated_row_size[i] ;
			for(size_t j = 0 ; j < A.row_size[i] ;j++)
			{
				if(A.column_index[start_index+j] != i)
					delta+=x[A.column_index[start_index+j]] * A.array[start_index+j] ;
			}
			double new_x_i = inverseDiagonal[i]*(b[i] - delta) ;
			
			
			err += std::abs(x[i]-new_x_i) ;
			x[i]=new_x_i ;
		}
		
// 		if(nit%1000 == 0 && verbose)
// 		{
// 		std::cout << "error :"<< err <<", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
// 		}
		
		nit++ ;
		
	}
// 	if(verbose)
// 		std::cout << "converged after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
	
	
	return x ;
}



}
