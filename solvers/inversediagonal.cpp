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

InverseDiagonal::InverseDiagonal(const CoordinateIndexedSparseMatrix &A) : diagonal(A.inverseDiagonal())
{
	double min = diagonal[0] ;
	double max = diagonal[0] ;
	for(size_t i = 1 ; i < diagonal.size() ; i++)
	{
		if(diagonal[i] < min)
			min = diagonal[i] ;

		if(diagonal[i] > max)
			max = diagonal[i] ;
	}
}

void  InverseDiagonal::precondition(const Vector &v, Vector & t) const
{
	for(size_t i = 0 ; i < t.size() ; i++)
		t[i]=v[i]*diagonal[i] ;
}

InverseDiagonalSquared::InverseDiagonalSquared(const CoordinateIndexedSparseMatrix &A)
{
	diagonal = new Vector(A.inverseDiagonalSquared()) ;
// 	double fac = std::abs(*diagonal).max()/std::abs(*diagonal).min();
// 	std::cout << "pseudo-C = " << fac << std::endl ;
}

void InverseDiagonalSquared::precondition(const Vector &v, Vector & t) const
{
	double * ti = &t[0] ;
	const double * vi = &v[0] ;
	#pragma omp parallel for
	for(int d = 0 ; d < diagonal->size() ; ++d)
		*(ti++)=(*vi++)*(*diagonal)[d] ;
}

