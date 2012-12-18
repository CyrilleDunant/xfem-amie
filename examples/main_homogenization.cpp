
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/matrixops.h"
#include "../physics/homogenization/composite.h"
#include "../physics/homogenization/phase.h"
#include "../physics/stiffness.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/materials/paste_behaviour.h"

#include <fstream>

using namespace Mu ;

Matrix i4_2D()
{
	Matrix ret(3,3) ;
	for(size_t i = 0 ; i < 3 ; i++)
	{
		ret[i][i] = 1. ;
		for(size_t j = 0 ; j < i ; j++)
		{
			ret[i][j] = 0. ;
			ret[j][i] = 0. ;
		}
	}
	return ret ;
}

Vector i2_2D()
{
	Vector ret(3) ;
	for(size_t i = 0 ; i < 2 ; i++)
	{
		ret[i] = 1. ;
	}
	return ret ;
}

Matrix invert(Matrix & m)
{
	return inverse3x3Matrix(m) ;
}

void print(Vector v)
{
	for(size_t i = 0 ; i < v.size() ; i++)
		std::cout << v[i] << "\t" ;
	std::cout << std::endl ;
}

int main(int argc, char *argv[])
{
	double vp = 0.3 ;
	double va = 0.7*0.99 ;
	double vg = 0.7*0.01 ;
	double fa = 11.e6 ;
	
	Matrix Cp = (new ElasticOnlyPasteBehaviour())->param ;
	Matrix Ca = (new ElasticOnlyAggregateBehaviour())->param ;
	Matrix Cg = (new GelBehaviour())->param ;
	
	Matrix Ap = i4_2D() ;
	Matrix Aa = i4_2D() ;
	Matrix Ag = i4_2D() ;
	
	Matrix Da = i4_2D() ;
	
	Vector bg = (Vector) (Cg*i2_2D())*(-0.22) ;
	
 	Matrix Ch = (Matrix) (Cp*Ap)*vp + (Matrix) (Da*Ca*Aa)*va + (Matrix) (Cg*Ag)*vg ;
 	Vector bh = (Vector) (Ag*bg)*vg ;
 	
	Vector sh = i2_2D()*0. ;
 	Vector eh = sh - (Vector) (invert(Ch)*bh) ;
	
	Vector ea = Aa*eh ;
	Vector sa = Da*Ca*ea ;
	bool met = false ;
	while(!met)
	{
		met = true ;
		for(size_t i = 0 ; i < 3 ; i++)
		{
			if(sa[i] > fa || sa[i] < fa*(-8))
			{
				met = false ;
				Da[i][i] -= 0.01 ;
			}
		}
		Ch = (Matrix) (Cp*Ap)*vp + (Matrix) (Da*Ca*Aa)*va + (Matrix) (Cg*Ag)*vg ;
		Ch.print() ;
		eh = sh - (Vector) (invert(Ch)*bh) ;
		ea = Aa*eh ;
		sa = (Da*Ca)*ea ;
//		print(sa) ;
	}
  
  
	return 0 ;
}
