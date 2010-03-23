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


Vector & LinearDamage::damageState() 
{
	return state ;
}
void LinearDamage::step(ElementState & s)
{
	double maxD = .99999 ;
	Vector pstrain = s.getPrincipalStresses(s.getParent()->getCenter()) ;
	Vector strain = s.getStrain(s.getParent()->getCenter()) ;
	
	bool inCompression = pstrain.min() < 0 && std::abs(pstrain.min()) > std::abs(pstrain.max()) ;

	double sum = 0 ;
	for(size_t i = 0 ; i < pstrain.size() ; i++)
		sum+=std::abs(pstrain[i]) ;
	double snorm = sqrt((pstrain*pstrain).sum()) ;
// 	std::cout << sum << " -> " << std::flush ;
	if(inCompression)
	{
		
// 		for(size_t i = 0 ; i < state.size()-1 ; i++)
// 		{
// 			state[i] += .1; //25e-5*(maxD/sqrt(s.getParent()->area()))*(std::abs(pstrain[i])/snorm) ; // ;
// 			state[i] = std::min(maxD, state[i]) ;
// 		}
		state[state.size()-1] += .1 ; //25e-5*(maxD/sqrt(s.getParent()->area())) ;
		state[state.size()-1] = std::min(maxD, state[state.size()-1]) ;
	}
	else/* if (!strainBroken)*/
	{
		for(size_t i = 0 ; i < state.size() ; i++)
		{
			state[i] += .1; //25e-5*(maxD/sqrt(s.getParent()->area()))*(std::abs(pstrain[i])/snorm) ; //5e-5*maxD/sqrt(s.getParent()->area()) ;
			state[i] = std::min(maxD, state[i]) ;
		}
	}
// 	std::cout << state.sum() << std::flush ;
}

void LinearDamage::artificialDamageStep(double d)
{
	for(size_t i = 0 ; i < state.size() -1 ; i++)
		state[i] = std::min(state[i]+d,0.9999) ;
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
	return state.max() >= .95 ;
}

LinearDamage::~LinearDamage()
{
}


}
