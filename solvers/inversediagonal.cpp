//
// C++ Implementation: inversediagonal
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "inversediagonal.h"

using namespace Mu ;


InverseLumpedDiagonal::InverseLumpedDiagonal(const CoordinateIndexedSparseMatrix &A)
{
	diagonal = new Vector(double(0), A.row_size.size()) ;
	
	size_t array_index = 0 ;
	
	for(size_t i = 0 ; i < A.row_size.size() ; i++)
	{
		for(size_t j =0 ; j < A.row_size[i] ; j++)
		{
			(*diagonal)[i] += A.array[array_index] ;
			array_index++ ;
		}
		
// 		if(std::abs((*diagonal)[i]) > 1e-12)
		(*diagonal)[i] = 1./(*diagonal)[i] ;
// 		else
// 			(*diagonal)[i] = 1000 ;
	}
}

void  InverseLumpedDiagonal::precondition(const Vector &v, Vector & t) const
{
	t=v*(*diagonal) ;
}

InverseDiagonal::InverseDiagonal(const CoordinateIndexedSparseMatrix &A)
{
	diagonal = new Vector(A.inverseDiagonal()) ;
// 	for(size_t i = 0 ; i < diagonal->size() ; i++)
// 	{
// 		if ((*diagonal)[i] < 0)
// 			(*diagonal)[i] = 1 ;
// 	}
// 	double fac = std::abs(*diagonal).max()/std::abs(*diagonal).min();
// 	std::cout << "pseudo-C = " << fac << std::endl ;
}

void  InverseDiagonal::precondition(const Vector &v, Vector & t) const
{
	t=v*(*diagonal) ;
}

InverseDiagonalSquared::InverseDiagonalSquared(const CoordinateIndexedSparseMatrix &A)
{
	diagonal = new Vector(A.inverseDiagonalSquared()) ;
// 	double fac = std::abs(*diagonal).max()/std::abs(*diagonal).min();
// 	std::cout << "pseudo-C = " << fac << std::endl ;
}

void InverseDiagonalSquared::precondition(const Vector &v, Vector & t) const
{
	t=v*(*diagonal) ;
}

