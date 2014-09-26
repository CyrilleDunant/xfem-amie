//
// C++ Implementation: choleskidecomposed
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "choleskidecomposed.h"

using namespace Amie ;

LowerTriangular::LowerTriangular(const CoordinateIndexedSparseMatrix &A_, const Vector &b_) :LinearSolver(A_, b_), d(A_.inverseDiagonal()) { };

bool LowerTriangular::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
    int stride = A.stride ;
	x = 0 ;
	
	for(size_t i = 0 ; i < A.row_size.size() ; i++)
    {
        for(size_t k = 0 ; k < stride ; k++)
        {
            double delta = 0 ; 
            int row = i*stride+k ;
            for(size_t j = 0 ; j < A.row_size[i] ; j++)
            {
                int column = A.column_index[A.accumulated_row_size[i]+j] ;
                for(size_t l = 0 ; l < stride ; l++)
                {
                    if(row >= column*stride+l)
                    {
                        delta += x[column*stride+l]*A[row][column*stride+l] ;
                    }
                }
                
                if(row < column*stride)
                    break ;
            }
            x[row] = (b[row] -delta)*d[row] ;
        }
    }
	
	return true ;
}

UpperTriangular::UpperTriangular(const CoordinateIndexedSparseMatrix &A_, const Vector &b_) :LinearSolver(A_, b_) , d(A_.inverseDiagonal()){ };

bool UpperTriangular::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
	x = 0 ;
    int stride = A.stride ;
	
	for(int i = A.row_size.size()-1 ; i >= 0 ; i--)
    {
        for(int k = stride-1 ; k >= 0 ; k--)
        {
            double delta = 0 ; 
            int row = i*stride+k ;
            for(int j = A.row_size[i]-1 ; j >= 0 ; j--)
            {
                int column = A.column_index[A.accumulated_row_size[i]+j] ;
                for(int l = stride-1 ; l >=0 ; l--)
                {
                    if(row <= column*stride+l)
                    {
                        delta += x[column*stride+l]*A[row][column*stride+l] ;
                    }
                }
                if(row > column*stride)
                    break ;
            }
            x[row] = (b[row] -delta)*d[row] ;
        }
    }
    
    return true ;

}

CholeskiDecomposed::CholeskiDecomposed(const CoordinateIndexedSparseMatrix &A_, const Vector &b_,const Vector &d_) :LinearSolver(A_, b_), d(d_), y(0., A_.row_size.size()*A_.stride)
{};

bool CholeskiDecomposed::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
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
	
	return true ;
}
