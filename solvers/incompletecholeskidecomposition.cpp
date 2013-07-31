
//
// C++ Implementation: incompletecholeskidecomposition
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "incompletecholeskidecomposition.h"
#include "choleskidecomposed.h"

using namespace Mu ;

InCompleteCholesky::InCompleteCholesky(const CoordinateIndexedSparseMatrix &A_) : A(A_), d(A_.diagonal())
{
	stable = true ;
	
	for(int i = 0 ; i < (int)A.row_size.size()*A.stride ; i++)
	{
		d[i] = 1./d[i] ;
		for(int j = i+1 ; j < (int)A.row_size.size()*A.stride ; j++)
			d[j] = A[i][j]*d[i] ;
		
		for(int j = i+1 ; j < (int)A.row_size.size()*A.stride ; j++)
		{
			for(int k = j ; k < (int)A.row_size.size()*A.stride ; k++)
			{
				if( std::abs(A[k][j]) < POINT_TOLERANCE_2D )
				{
					double dummy = -A[i][j]*d[k] ;
					A[k][k] += std::abs(dummy) ;
					A[j][j] += std::abs(dummy) ;
				}
				else
				{
					double dummy = A[i][j] ;
					A[k][j] -= d[k]*dummy ;
					A[j][k] -= d[k]*dummy ;
				}
			}
		}
	}
}

void InCompleteCholesky::precondition(const Vector& v, Vector& t)
{
// 	if(stable)
// 	{
//		Vector y = UpperTriangular(A, v).solve(Vector(0), nullptr) ;
//		Vector ret= LowerTriangular(A, y).solve(Vector(0), nullptr) ;
	Vector v_(v) ;
	CholeskiDecomposed sv(A,v_,d) ;
	sv.solve(Vector(0), nullptr) ;
	t = sv.x ;
// 	}
// 	else
// 		return UpperTriangular(A, v).solve(Vector(0), nullptr) ;
}
