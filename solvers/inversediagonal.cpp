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

#include "inversediagonal.h"
#include <limits>
#include <iostream>

using namespace Amie ;

InverseLumpedDiagonal::InverseLumpedDiagonal(const CoordinateIndexedSparseMatrix &A) : diagonal(0., A.row_size.size()*A.stride)
{

	for(size_t i = 0 ; i < A.row_size.size()*A.stride ; i++)
	{

		diagonal[i] += A[i][0] ;

		
		for(size_t j = 1 ; j < A.row_size.size()*A.stride ; j++)
		{
			diagonal[i] += A[j][i] ;
		}
		
		double v = diagonal[i] ;
		if(std::abs(v) > std::numeric_limits<double>::epsilon())
			diagonal[i] = 1./v ;
		else if(v > 0)
			diagonal[i] = 1./std::numeric_limits<double>::epsilon() ;
        else
            diagonal[i] = -1./std::numeric_limits<double>::epsilon() ;

		
// 		std::cout << diagonal[i] << std::endl ;
	}
}

void  InverseLumpedDiagonal::precondition(const Vector &v, Vector & t) 
{
	for(size_t i = 0 ; i < v.size() ; i++)
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
	#pragma omp parallel for schedule(static) if (t.size() > 10000)
	for(size_t i = 0 ; i < v.size() ; i++)
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
	#pragma omp parallel for schedule(static) if (t.size() > 10000)
	for(size_t d = 0 ; d < diagonal->size() ; ++d)
		*(ti++)=(*vi++)*(*diagonal)[d] ;
}

Inverse2x2Diagonal::Inverse2x2Diagonal(const CoordinateIndexedSparseMatrix &A)
{
	if(A.row_size.size()*A.stride%2 == 1)
	{
		std::cerr << "unable to create Inverse2x2Diagonal preconditionner with odd number of dof" << std::endl ;
		exit(0) ;
	}

	for(size_t i = 0 ; i < A.row_size.size()*A.stride ; i += 2)
	{
		Matrix block(2,2) ;
		block[0][0] = A[i][i] ;
		block[1][0] = A[i+1][i] ;
		block[0][1] = A[i][i+1] ;
		block[1][1] = A[i+1][i+1] ;
		if(std::abs(det(block)) > std::numeric_limits<double>::epsilon())
		{
			invert2x2Matrix(block) ;
		}
		else
		{
			block[0][1] = 0. ;
			block[1][0] = 0. ;
			if(std::abs(A[i][i]) > std::numeric_limits<double>::epsilon())
				block[0][0] = 1./A[i][i] ;
			else
				block[0][0] = 1. ;
			if(std::abs(A[i+1][i+1]) > std::numeric_limits<double>::epsilon())
				block[1][1] = 1./A[i+1][i+1] ;
			else
				block[1][1] = 1. ;
		}
		blocks.push_back(block) ;
	}

}

void  Inverse2x2Diagonal::precondition(const Vector &v, Vector & t) 
{
	if(blocks.size()*2 != v.size())
	{
		std::cout << "preconditionner and vector do not have the same size!" << std::endl ;
	}

	for(size_t i = 0 ; i < blocks.size() ; i++)
	{
		t[i*2] = blocks[i][0][0]*v[i*2]+blocks[i][0][1]*v[i*2+1] ;
		t[i*2+1] = blocks[i][1][0]*v[i*2]+blocks[i][1][1]*v[i*2+1];
	}
}

