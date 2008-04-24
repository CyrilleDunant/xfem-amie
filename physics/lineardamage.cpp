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
	Vector pstrain = s.getPrincipalStresses(s.getParent()->getCenter()) ;
	Vector strain = s.getStrain(s.getParent()->getCenter()) ;
	
	bool inCompression = pstrain.min() < 0 && std::abs(pstrain.min()) > pstrain.max() ;
	bool strainBroken = strain.max() > strainLimit ;

	double sum = 0 ;
	for(size_t i = 0 ; i < pstrain.size() ; i++)
		sum+=std::abs(pstrain[i]) ;

	if(inCompression)
	{
		std::cout << "crunch..." << std::endl ;
		state[state.size()-1] += std::max(.01, state[state.size()-1]) ;
		state[state.size()-1] = std::min(.9999, state[state.size()-1]) ;
	}
	else if (!strainBroken)
	{
		std::cout << "crack..." << std::endl ;
		for(size_t i = 0 ; i < state.size()-1 ; i++)
		{
			if(sum > 1e-12)
				state[i] += std::max(.01, state[i])*std::abs(pstrain[i])/sum ;
			else
				state[i] += .9 ;
			state[i] = std::min(.9999, state[i]) ;
		}
	}
	else if (strainBroken)
	{
		std::cout << "crack!" << std::endl ;
		for(size_t i = 0 ; i < state.size() ; i++)
			state[i] = .9999 ;
	}
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
