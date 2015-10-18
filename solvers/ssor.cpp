//
// C++ Implementation: SSOR preconditinonner
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


Ssor::Ssor(const CoordinateIndexedSparseMatrix &A, int rowstart , int colstart ) : upper(A), lower(A), buffer(A.inverseDiagonal()), omega(1.5), rowstart(rowstart), colstart(colstart)
{
     
    //first, update the sparseness patterns

    int stride = A.stride ;
    std::set<std::pair<size_t, size_t>> upperPattern ;
    std::set<std::pair<size_t, size_t>> lowerPattern ;
    size_t column_iterator = 0 ;
    for(size_t i = 0 ; i < A.row_size.size() ; i++)
    {
        size_t row = i ;
        for(size_t j = 0 ; j < A.row_size[i] ; j++)
        {
            
            size_t column = A.column_index[column_iterator++] ;
            
            if (row <= column)
                upperPattern.insert(std::make_pair(row,column)) ;
            if (row >= column)
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
        for(int k = 0 ; k < stride ; k++)
        {                   
            double v = A[row * stride+k][row * stride+k] ;
            long double diag = 0 ;         
            
            if (v > 1e-12)
            {
                diag = 1./v ;
            }
            
            for(size_t j = 0 ; j < A.row_size[i] ; j++)
            {
                int column = A.column_index[A.accumulated_row_size[i]+j] ;

                for(int l = 0 ; l < stride ; l++)
                {    
                    if(row*stride+k > column*stride+l)
                    {
                        long double val = -sqrt((2.-omega)*omega)*sqrt(diag)*A[row*stride+k][column*stride+l]*diag*omega;
                        upper[row*stride+k][column*stride+l] = val ;
                        lower[column*stride+l][row*stride+k] = val ;
                        
                    }
                    else if(row*stride+k == column*stride+l)
                    {
                        long double val = sqrt((2.-omega)*omega)*sqrt(diag);
                        upper[row*stride+k][column*stride+l] = val ;
                        lower[column*stride+l][row*stride+k] = val ;
                    }
                }
            }
        }
    }
    
}
// y = (upper+auxupper*auxuper)*(lower+auxlower*auxlower)*b
//   = (upper+auxupper*auxupper)*(lower*b+auxlower*auxlower*b)
//   =  upper*lower*b + upper*auxlower*auxlower*b + auxupp*auxupper*lower*b + auxupper*auxupper*auxlower*auxlower*b
void  Ssor::precondition(const Vector &v, Vector & t) 
{
    if(v.size() != buffer.size())
      std::cout << " ouch ! " << std::endl ;

    assign(buffer, upper*v, rowstart, colstart) ;
    assign(t, lower*buffer, rowstart, colstart) ; //1

    
}

