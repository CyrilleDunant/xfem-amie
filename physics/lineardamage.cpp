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
 : state(numDof + 1),strainLimit(threshold)
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
	double maxD = .9999 ;
	Vector pstrain = s.getPrincipalStresses(s.getParent()->getCenter()) ;
	Vector strain = s.getStrain(s.getParent()->getCenter()) ;
	
	bool inCompression = pstrain.min() < 0 && std::abs(pstrain.min()) > pstrain.max() ;
	bool strainBroken = strain.max() > strainLimit ;

	double sum = 0 ;
	for(size_t i = 0 ; i < pstrain.size() ; i++)
		sum+=std::abs(pstrain[i]) ;

	if(inCompression)
	{
		state[state.size()-1] += 5e-5*maxD/sqrt(s.getParent()->area()) ;
		state[state.size()-1] = .5 ;
	}
	else if (!strainBroken)
	{
		for(size_t i = 0 ; i < state.size()-1 ; i++)
		{
			if(sum > 1e-12)
				state[i] += 5e-5*maxD/sqrt(s.getParent()->area())*std::abs(pstrain[i])/sum ;
			else
				state[i] += 5e-5*maxD/sqrt(s.getParent()->area()) ;
			state[i] = std::min(maxD, state[i]) ;
		}
	}
	else if (strainBroken)
	{
		for(size_t i = 0 ; i < state.size() ; i++)
			state[i] = maxD ;
	}
}

Matrix LinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;
	
	for(size_t i = 0 ; i < (m.numRows()+1)/2 ;i++)
	{
		for(size_t j = 0 ; j < m.numCols() ;j++)
		{
			ret[i][j] *= 1.-(state[0]) ;
		}
	}

	for(size_t i = (m.numRows()+1)/2 ; i < m.numRows() ;i++)
	{
		for(size_t j = 0 ; j < m.numCols() ;j++)
		{
			ret[i][j] *= 1.-(state[i]) ;
		}
	}

	return ret ;
}


bool LinearDamage::fractured() const
{
	return state.max() >= .999999 ;
}

LinearDamage::~LinearDamage()
{
}


}
