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
	diagonal = new Vector(0., A.row_size.size()*(A.stride+A.stride%2)) ;
	Vector eye(1., A.row_size.size()*(A.stride+A.stride%2)) ;
	
	size_t array_index = 0 ;

	for(size_t i = 0 ; i < diagonal->size() ; i += (A.stride+A.stride%2) )
	{
		Vector v = A[i]*eye ;
		for(size_t j = 0 ; j < (A.stride+A.stride%2) ; j++)
			(*diagonal)[i+j] = 1./v[j] ;
	}
}

void  InverseLumpedDiagonal::precondition(const Vector &v, Vector & t) 
{
	t=v*(*diagonal) ;
}

InverseDiagonal::InverseDiagonal(const CoordinateIndexedSparseMatrix &A) : diagonal(A.inverseDiagonal())
{
// 	double min = diagonal[0] ;
// 	double max = diagonal[0] ;
// 	for(size_t i = 1 ; i < diagonal.size() ; i++)
// 	{
// 		if(diagonal[i] < min)
// 			min = diagonal[i] ;
// 
// 		if(diagonal[i] > max)
// 			max = diagonal[i] ;
// 	}
}

void  InverseDiagonal::precondition(const Vector &v, Vector & t) 
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

void InverseDiagonalSquared::precondition(const Vector &v, Vector & t) 
{
	double * ti = &t[0] ;
	const double * vi = &v[0] ;
	#pragma omp parallel for
	for(int d = 0 ; d < diagonal->size() ; ++d)
		*(ti++)=(*vi++)*(*diagonal)[d] ;
}

