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


Ssor::Ssor(const CoordinateIndexedSparseMatrix &A, double omega) : upper(A), lower(A),  buffer(1., A.row_size.size()*A.stride), omega(omega), omega_prev(0), previous_r(0)
{
     
    //first, update the sparseness patterns
    size_t blocksUpper = 0 ;
    size_t blocksLower = 0 ;
    int stride = A.stride ;
    for(size_t i = 0 ; i < A.row_size.size() ; i++)
    {
        for(size_t j = 0 ; j < A.row_size[i] ; j++)
        {
            int row = i ;
            int column = A.column_index[A.accumulated_row_size[i]+j] ;
            
            if(row > column)
                blocksLower++ ;
            else if (row < column)
                blocksUpper++ ;
            else
            {
                blocksLower++ ;
                blocksUpper++ ;
            }    
        }
    }
    upper.column_index.resize(blocksUpper) ;
    lower.column_index.resize(blocksLower) ;
    upper.array.resize(blocksUpper*stride*(stride+stride%2), 0.) ; 
    lower.array.resize(blocksUpper*stride*(stride+stride%2), 0.) ; 
    
    int upperPreviousRowSize = 0;
    int lowerPreviousRowSize = 0;
    int lowerIterator = 0 ;
    int upperIterator = 0;
    for(size_t i = 0 ; i < A.row_size.size() ; i++)
    {
        if(i > 0)
        {
            upper.accumulated_row_size[i] = upperPreviousRowSize+upper.accumulated_row_size[i-1] ;
            lower.accumulated_row_size[i] = lowerPreviousRowSize+lower.accumulated_row_size[i-1] ;
        }
        upperPreviousRowSize = 0;
        lowerPreviousRowSize = 0;
        for(size_t j = 0 ; j < A.row_size[i] ; j++)
        {
            int row = i ;
            int column = A.column_index[A.accumulated_row_size[i]+j] ;
            
            if(row > column)
            {
                lower.column_index[lowerIterator++] = column ;
                lowerPreviousRowSize++ ;
            }
            else if (row < column)
            {
                upper.column_index[upperIterator++] = column ;
                upperPreviousRowSize++ ;
            }
            else
            {
                lower.column_index[lowerIterator++] = column ;
                lowerPreviousRowSize++ ;
                upper.column_index[upperIterator++] = column ;
                upperPreviousRowSize++ ;
            }    
        }
        upper.row_size[i] = upperPreviousRowSize ;
        lower.row_size[i] = lowerPreviousRowSize ;
    }

    
    for(size_t i = 0 ; i < A.row_size.size() ; i++)
    {
        for(size_t j = 0 ; j < A.row_size[i] ; j++)
        {
            int row = i ;
            int column = A.column_index[A.accumulated_row_size[i]+j] ;
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
                        upper[row*stride+k][column*stride+l] = val ;
                        lower[column*stride+l][row*stride+k] = val ;
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

    assign(buffer, lower*v) ;
    assign(t, upper*buffer) ;
    
}

