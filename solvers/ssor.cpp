//
// C++ Implementation: inversediagonal
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "ssor.h"
#include "choleskidecomposed.h"
#include <limits>
#include <iostream>

using namespace Amie ;


Ssor::Ssor(const CoordinateIndexedSparseMatrix &A, double omega) : upper(A), lower(A),  buffer(A.inverseDiagonal()), omega(omega), omega_prev(0), previous_r(0)
{
     
    //first, update the sparseness patterns
    size_t blocksUpper = 0 ;
    size_t blocksLower = 0 ;
    int stride = A.stride ;
    std::set<std::pair<size_t, size_t>> upperPattern ;
    std::set<std::pair<size_t, size_t>> lowerPattern ;
    for(size_t i = 0 ; i < A.row_size.size() ; i++)
    {
        int row = i ;
        for(size_t j = 0 ; j < A.row_size[i] ; j++)
        {
            
            int column = A.column_index[A.accumulated_row_size[i]+j] ;
            
            if (row >= column)
                upperPattern.insert(std::make_pair(row,column)) ;
            if (row <= column)
                lowerPattern.insert(std::make_pair(row,column)) ;
        }
    }
    for(size_t i = 0 ; i < A.row_size.size() ; i++)
    {
       upperPattern.insert(std::make_pair(i,i)) ;
       lowerPattern.insert(std::make_pair(i,i)) ;
    }
    upper.reshape(upperPattern, stride) ;
    lower.reshape(lowerPattern, stride) ;

    for(size_t i = 0 ; i < A.row_size.size() ; i++)
    {
        int row = i ;
        for(size_t j = 0 ; j < A.row_size[i] ; j++)
        {
            
            int column = A.column_index[A.accumulated_row_size[i]+j] ;
	    if(row < column)
	      continue ;
            for(size_t k = 0 ; k < stride ; k++)
            {
                double v = A[row * stride+k][row * stride+k] ;
                if (std::abs(v) < 1e-8)
                    continue ;
                double diag = 1./v ;

                for(size_t l = 0 ; l < stride ; l++)
                {
                     
                    if(row*stride+k > column*stride+l)
                    {
                        double val = -sqrt(diag)*A[row*stride+k][column*stride+l]*diag;
// 			std::cout << "a" << std::endl ;
                        upper[row*stride+k][column*stride+l] = val ;
// 			std::cout << "b" << std::endl ;
                        lower[column*stride+l][row*stride+k] = val ;
// 			std::cout << "b" << std::endl ;
			
                    }
                    else if(row*stride+k == column*stride+l)
                    {
                        double val = sqrt(diag);
                        upper[row*stride+k][column*stride+l] = val ;
                        lower[column*stride+l][row*stride+k] = val ;
                    }
                }
            }
        }
    }
    upper.array *= sqrt((2.-omega)*omega)*omega ;
    lower.array *= sqrt((2.-omega)*omega)*omega ;
    
    for(size_t i = 0 ; i < A.row_size.size() ; i++)
    {
        for(size_t k = 0 ; k < stride ; k++)
        {
            upper[i*stride+k][i*stride+k] /= omega ;
            lower[i*stride+k][i*stride+k] /= omega ;
        }
    }
}

void  Ssor::precondition(const Vector &v, Vector & t) 
{
    if(v.size() != buffer.size())
      std::cout << " ouch ! " << std::endl ;
    assign(buffer, lower*v) ;
    assign(t, upper*buffer) ;
    
}

