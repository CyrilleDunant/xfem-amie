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
#include <limits>
#include <iostream>

using namespace Mu ;

InverseLumpedDiagonal::InverseLumpedDiagonal(const CoordinateIndexedSparseMatrix &A) : diagonal(0., A.row_size.size()*A.stride)
{

	for(size_t i = 0 ; i < A.row_size.size()*A.stride ; i++)
	{
		for(size_t j = 0 ; j < 1 ; j++)
		{
			diagonal[i] += A[i][j] ;
		}
		
		for(size_t j = 1 ; j < A.row_size.size()*A.stride ; j++)
		{
			diagonal[i] += A[j][i] ;
		}
		
		double v = diagonal[i] ;
		if(std::abs(v) > std::numeric_limits<double>::epsilon())
			diagonal[i] = 1./v ;
		else
			diagonal[i] = 1./std::numeric_limits<double>::epsilon() ;

		
		std::cout << diagonal[i] << std::endl ;
	}
}

void  InverseLumpedDiagonal::precondition(const Vector &v, Vector & t) 
{
	for(size_t i = 0 ; i < t.size() ; i++)
		t[i]=v[i]*diagonal[i] ;
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

