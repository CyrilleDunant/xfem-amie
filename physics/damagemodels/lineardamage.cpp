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

LinearDamage::LinearDamage(int numDof)
 : state(numDof + 1)
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
	
	bool inCompression = pstrain.min() < 0 && std::abs(pstrain.min()) > std::abs(pstrain.max()) ;

	double sum = 0 ;
	for(size_t i = 0 ; i < pstrain.size() ; i++)
		sum+=std::abs(pstrain[i]) ;

	if(inCompression)
	{
		double snorm = sqrt((pstrain*pstrain).sum()) ;
		for(size_t i = 0 ; i < state.size()-1 ; i++)
		{
			state[i] += .005*std::abs(pstrain[i])/snorm ; //5e-5*maxD/sqrt(s.getParent()->area()) ;
			state[i] = std::min(maxD, state[i]) ;
		}
		state[state.size()-1] += .1 ;
		state[state.size()-1] = std::min(maxD, state[state.size()-1]) ;
	}
	else/* if (!strainBroken)*/
	{
		for(size_t i = 0 ; i < state.size()-1 ; i++)
		{
			state[i] += .01*std::abs(pstrain[i])/std::abs(pstrain).sum() ; //5e-5*maxD/sqrt(s.getParent()->area()) ;
			state[i] = std::min(maxD, state[i]) ;
		}
	}

}

Matrix LinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;
	
	if(fractured())
		return m*0.0001	;
	//this is a silly way of distinguishing between 2D and 3D
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
	return state.max() >= .85 ;
}

LinearDamage::~LinearDamage()
{
}


}
