//
// C++ Implementation: choleskidecomposed
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "choleskidecomposed.h"

using namespace Mu ;

LowerTriangular::LowerTriangular(const CoordinateIndexedSparseMatrix &A_, const Vector &b_) :LinearSolver(A_, b_), d(A_.inverseDiagonal()) { };

Vector & LowerTriangular::solve(const Vector &x0, const Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
	x = 0 ;
	for(size_t i = 0 ; i < x.size() ; i++)
	{
		double delta  = 0 ;
		for(size_t j = A.accumulated_row_size[i] ; A.column_index[j]  < i ; j++)
		{
			delta += x[A.column_index[j]]*A.array[j] ;
		}
		x[i] = (b[i] -delta)*d[i] ;
	}
	
	return x ;
}

UpperTriangular::UpperTriangular(const CoordinateIndexedSparseMatrix &A_, const Vector &b_) :LinearSolver(A_, b_) , d(A_.inverseDiagonal()){ };

Vector & UpperTriangular::solve(const Vector &x0, const Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
	x = 0 ;
	for(size_t i = x.size() ; i > 0 ; i--)
	{
		double delta  = 0 ;
		for(size_t j = A.accumulated_row_size[i-1]+A.row_size[i-1] ; A.column_index[j-1]  > i-1 ; j--)
		{
			delta += x[A.column_index[j-1]]*A.array[j-1] ;
		}
		x[i-1] = (b[i-1] -delta)*d[i-1] ;
	}
	
	return x ;
}

CholeskiDecomposed::CholeskiDecomposed(const CoordinateIndexedSparseMatrix &A_, const Vector &b_,const Vector &d_) :LinearSolver(A_, b_), d(d_), y(0., A_.row_size.size())
{};

Vector & CholeskiDecomposed::solve(const Vector &x0, const Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
	x = 0 ;
	y = 0 ;
	for(int i = x.size()-1 ; i >= 0 ; --i)
	{
// 		x[i] = ( b[i] - reverseInnerProduct(A[i], x, i))*d[i] ;
		reverseInnerProductAssignAndAdd(A[i], x,x[i], b[i] ,i) ;
		x[i]*=-d[i] ;
	}
	for(size_t i = 0 ; i < x.size() ; ++i)
	{
// 		y[i] = (x[i] -innerProduct(A[i], y, i))*d[i];
		innerProductAssignAndAdd(A[i], y,y[i], x[i] ,i) ;
		y[i]*=-d[i] ;
	}
	
	return y ;
}
