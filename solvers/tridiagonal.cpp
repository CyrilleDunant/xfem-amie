//
// C++ Implementation: inversediagonal
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "tridiagonal.h"

using namespace Mu ;


TriDiagonal::TriDiagonal(const CoordinateIndexedSparseMatrix &A) : diagonal(A.diagonal()), upper(diagonal.size()-1)
{
	for(size_t i = 0 ; i < upper.size() ; i++)
		upper[i] = A[i][i+1] ;
	for(size_t i = 0 ; i < diagonal.size() ; i++)
	{
		if(std::abs(diagonal[i]) < 1e-12)
		{
			diagonal[i] = 1e-12 ;
			std::cout << "\n!!" << std::endl ;
		}
	}
}

void  TriDiagonal::precondition(const Vector &force, Vector & solution)
{
	Vector b(diagonal) ;
	Vector c(upper) ;
	Vector d(force) ;
	int n = diagonal.size() ;
	c[0] /= b[0];
	d[0] /= b[0];
	for(int i = 1 ; i < n-1 ; i++)
	{
		double id = 1.0/(b[i] - c[i - 1]*upper[i-1]); 
		if(i < n-1)
			c[i] *= id;
		d[i] = (d[i] - d[i - 1]*upper[i-1])*id;
	}
	
	d[n-1] = (d[n-1] - d[n - 2]*upper[n-2])/(b[n-1] - c[n - 2]*upper[n-2]);
	
	solution[n - 1] = d[n - 1];
	for(int i = n - 2; i >= 0; i--)
		solution[i] = d[i] - c[i]*solution[i + 1];
}

