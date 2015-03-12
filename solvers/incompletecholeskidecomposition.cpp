
//
// C++ Implementation: incompletecholeskidecomposition
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "incompletecholeskidecomposition.h"
#include "choleskidecomposed.h"

using namespace Amie ;

InCompleteCholesky::InCompleteCholesky(Assembly * A_) : d(A_->getMatrix().diagonal()), A(A_)
{
	stable = true ;
	
	for(int i = 0 ; i < (int)A->getMatrix().row_size.size()*A->getMatrix().stride ; i++)
	{
		d[i] = 1./d[i] ;
		for(int j = i+1 ; j < (int)A->getMatrix().row_size.size()*A->getMatrix().stride ; j++)
			d[j] = A->getMatrix()[i][j]*d[i] ;
		
		for(int j = i+1 ; j < (int)A->getMatrix().row_size.size()*A->getMatrix().stride ; j++)
		{
			for(int k = j ; k < (int)A->getMatrix().row_size.size()*A->getMatrix().stride ; k++)
			{
				if( std::abs(A->getMatrix()[k][j]) < POINT_TOLERANCE_2D )
				{
					double dummy = -A->getMatrix()[i][j]*d[k] ;
					A->getMatrix()[k][k] += std::abs(dummy) ;
					A->getMatrix()[j][j] += std::abs(dummy) ;
				}
				else
				{
					double dummy = A->getMatrix()[i][j] ;
					A->getMatrix()[k][j] -= d[k]*dummy ;
					A->getMatrix()[j][k] -= d[k]*dummy ;
				}
			}
		}
	}
}

void InCompleteCholesky::precondition(const Vector& v, Vector& t)
{
// 	if(stable)
// 	{
//		Vector y = UpperTriangular(A->getMatrix(), v).solve(Vector(0), nullptr) ;
//		Vector ret= LowerTriangular(A->getMatrix(), y).solve(Vector(0), nullptr) ;
	Vector v_(v) ;
	CholeskiDecomposed sv(A, v_) ;
	sv.solve(Vector(0), nullptr) ;
	t = sv.x ;
// 	}
// 	else
// 		return UpperTriangular(A->getMatrix(), v).solve(Vector(0), nullptr) ;
}
