//
// C++ Implementation: lineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "lineardamage.h"

namespace Mu {

LinearDamage::LinearDamage(int numDof, double threshold)
 : state(numDof + 1),thresholdDensity(threshold)
{
	isNull = false ;
	state = 0 ;
}

const Vector & LinearDamage::damageState() const
{
	return state ;
}

void LinearDamage::step(ElementState & s)
{
	Vector pstrain = s.getPrincipalStresses(s.getParent()->getCenter()) ;
	
	bool inCompression = false ;
	for(size_t i = 0 ; i < state.size()-1 ; i++)
	{
		if(pstrain[i] < 0)
		{
			inCompression = true ;
			break ;
		}
	}

	double sum = 0 ;
	for(size_t i = 0 ; i < pstrain.size() ; i++)
		sum+=std::abs(pstrain[i]) ;

	if(inCompression)
	{
		state[state.size()-1] += .15 ;
		state[state.size()-1] = std::min(.999, state[state.size()-1]) ;
	}
// 	else
// 	{
		for(size_t i = 0 ; i < state.size()-1 ; i++)
		{
			if(sum > 1e-12)
				state[i] += .05*std::abs(pstrain[i])/sum ;
			else
				state[i] += .05 ;
			state[i] = std::min(.999, state[i]) ;
		}
// 	}
}

Matrix LinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;
	
	for(size_t i = 0 ; i < m.numRows() ;i++)
	{
		for(size_t j = 0 ; j < m.numCols() ;j++)
		{
			ret[i][j] *= 1.-sqrt((state[i])*(state[j])) ;
		}
	}

	return ret ;
}


LinearDamage::~LinearDamage()
{
}


}
