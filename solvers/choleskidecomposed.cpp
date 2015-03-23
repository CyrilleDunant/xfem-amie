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

LowerTriangular::LowerTriangular(Assembly * a) :LinearSolver(a), d(a->getMatrix().inverseDiagonal()) { }

bool LowerTriangular::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
    int stride = assembly->getMatrix().stride ;
    x = 0 ;

    for(size_t i = 0 ; i < assembly->getMatrix().row_size.size() ; i++)
    {
        for(int k = 0 ; k < stride ; k++)
        {
            double delta = 0 ;
            int row = i*stride+k ;
            for(size_t j = 0 ; j < assembly->getMatrix().row_size[i] ; j++)
            {
                int column = assembly->getMatrix().column_index[assembly->getMatrix().accumulated_row_size[i]+j] ;
                for(int l = 0 ; l < stride ; l++)
                {
                    if(row >= column*stride+l)
                    {
                        delta += x[column*stride+l]*assembly->getMatrix()[row][column*stride+l] ;
                    }
                }

                if(row < column*stride)
                    break ;
            }
            x[row] = (assembly->getForces()[row] -delta)*d[row] ;
        }
    }

    return true ;
}

UpperTriangular::UpperTriangular(Assembly * a) :LinearSolver(a) , d(a->getMatrix().inverseDiagonal()) { }

bool UpperTriangular::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
    x = 0 ;
    int stride = assembly->getMatrix().stride ;

    for(int i = assembly->getMatrix().row_size.size()-1 ; i >= 0 ; i--)
    {
        for(int k = stride-1 ; k >= 0 ; k--)
        {
            double delta = 0 ;
            int row = i*stride+k ;
            for(int j = assembly->getMatrix().row_size[i]-1 ; j >= 0 ; j--)
            {
                int column = assembly->getMatrix().column_index[assembly->getMatrix().accumulated_row_size[i]+j] ;
                for(int l = stride-1 ; l >=0 ; l--)
                {
                    if(row <= column*stride+l)
                    {
                        delta += x[column*stride+l]*assembly->getMatrix()[row][column*stride+l] ;
                    }
                }
                if(row > column*stride)
                    break ;
            }
            x[row] = (assembly->getForces()[row] -delta)*d[row] ;
        }
    }

    return true ;

}

CholeskiDecomposed::CholeskiDecomposed(Assembly *a,const Vector &d_) :LinearSolver(a), d(d_), y(0., a->getMatrix().row_size.size()*a->getMatrix().stride)
{}

bool CholeskiDecomposed::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
    x = 0 ;
    y = 0 ;
    for(int i = x.size()-1 ; i >= 0 ; --i)
    {
// 		x[i] = ( b[i] - reverseInnerProduct(A[i], x, i))*d[i] ;
        reverseInnerProductAssignAndAdd(assembly->getMatrix()[i], x,x[i], assembly->getForces()[i] ,i) ;
        x[i]*=-d[i] ;
    }
    for(size_t i = 0 ; i < x.size() ; ++i)
    {
// 		y[i] = (x[i] -innerProduct(A[i], y, i))*d[i];
        innerProductAssignAndAdd(assembly->getMatrix()[i], y,y[i], x[i] ,i) ;
        y[i]*=-d[i] ;
    }

    return true ;
}
