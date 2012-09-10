//
// C++ Implementation: lineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "anisotropicdamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu {

AnisotropicLinearDamage::AnisotropicLinearDamage() 
{
// 	getState(true).resize(3, 0.) ;
// 	getPreviousState().resize(3, 0.) ;
	isNull = false ;

}

std::pair<Vector, Vector> AnisotropicLinearDamage::computeDamageIncrement(ElementState & s)
{

	
	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL && getState().size() == 0)
	{
		upState.resize(2, 0.) ;
		downState.resize(2, 0.) ;
		getState(true).resize(2, 0.) ;
	}
	else if( getState().size() == 0 && s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	{
		upState.resize(3, 0.) ;
		downState.resize(3, 0.) ;
		getState(true).resize(3, 0.) ;
	}
	
	Vector ret(0., state.size()) ;
// 	double E_2 = s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())[0][0] ; E_2*=E_2 ;
// 	double l_2 = s.getParent()->area() ; 
// 	double maxincrement = std::abs((l_2*E_2-1.)/(l_2+l_2*E_2)) ;
	double tensionDamagex = 0;
	double tensionDamagey = 0;
	double tensionDamagez = 0;
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression(1))
	{
		ret = 1 ;
	}
	else 
	{
		Vector stress(0., 3+3*(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)) ;
		s.getAverageField( REAL_STRESS_FIELD, stress) ;
		
		if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		{
			tensionDamagex = stress[0]*(stress[0]>0) ;//std::abs(cos(angle[0])) ;
			tensionDamagey = stress[1]*(stress[1]>0) ;//std::abs(cos(angle[1])) ;
			tensionDamagez = stress[2]*(stress[2]>0) ;//std::abs(cos(angle[2])) ;
			ret[2] = tensionDamagez*(getState()[2] < thresholdDamageDensity);
		}
		else
		{
			tensionDamagex = stress[0]*(stress[0]>0) ;//std::abs(cos(angle[0])) ;
			tensionDamagey = stress[1]*(stress[1]>0) ;//std::abs(sin(angle[0])) ;
		}
		ret[0] = tensionDamagex*(getState()[0] < thresholdDamageDensity);
		ret[1] = tensionDamagey*(getState()[1] < thresholdDamageDensity);
		
	}
	if(ret.max() > POINT_TOLERANCE_2D)
		ret /= ret.max() ;
	Vector factor = -(state-1.) ;
	
	for(int i = 0 ; i < ret.size() ; i++)
		if(ret[i] < POINT_TOLERANCE_2D)
			ret[i] = 1 ;
		
	factor /= ret ;
	return std::make_pair(state,state+ret*factor.min()) ;
// 	std::cout << getState().sum() << std::flush ;
}

void AnisotropicLinearDamage::computeDelta(const ElementState & s)
{

	
	Vector ret = state ;
	// 	double E_2 = s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())[0][0] ; E_2*=E_2 ;
	// 	double l_2 = s.getParent()->area() ; 
	// 	double maxincrement = std::abs((l_2*E_2-1.)/(l_2+l_2*E_2)) ;
	double tensionDamagex = 0;
	double tensionDamagey = 0;
	double tensionDamagez = 0;
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression(1))
	{
		ret = 1 ;
	}
	else 
	{
		Vector stress(0., 3+3*(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)) ;
		Point center = s.getParent()->getCenter() ;
		s.getField( REAL_STRESS_FIELD, center, stress, false) ;
		
		if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		{
			tensionDamagex = stress[0]*(stress[0]>0) ;//std::abs(cos(angle[0])) ;
			tensionDamagey = stress[1]*(stress[1]>0) ;//std::abs(cos(angle[1])) ;
			tensionDamagez = stress[2]*(stress[2]>0) ;//std::abs(cos(angle[2])) ;
			ret[2] = tensionDamagez*(getState()[2] < thresholdDamageDensity);
		}
		else
		{
			tensionDamagex = stress[0]*(stress[0]>0) ;//std::abs(cos(angle[0])) ;
			tensionDamagey = stress[1]*(stress[1]>0) ;//std::abs(sin(angle[0])) ;
		}
		ret[0] = tensionDamagex*(getState()[0] < thresholdDamageDensity);
		ret[1] = tensionDamagey*(getState()[1] < thresholdDamageDensity);
		
	}
	if(ret.max() > POINT_TOLERANCE_2D)
		ret /= ret.max() ;
	Vector factor = -(state-1.) ;
	
	for(int i = 0 ; i < ret.size() ; i++)
		if(ret[i] < POINT_TOLERANCE_2D)
			ret[i] = 1 ;
		
		factor /= ret ;
	delta = (state+ret*factor.min()-state).max() ;
}


Matrix AnisotropicLinearDamage::apply(const Matrix & m, const Point & p, const IntegrableEntity * e, int g ) const
{
	if(state.max() < POINT_TOLERANCE_2D)
		return m ;
	Matrix ret(m) ;
	
	if(fractured())
		return m*0.;
	Vector rstate = getState() ;
	for(size_t i = 0 ; i < getState().size() ; i++)
	{
		rstate[i] = std::min(thresholdDamageDensity,rstate[i]) ;
	}

	ret[0][0] = m[0][0]*(1.-rstate[0]) ; ret[0][1] = m[0][1]*sqrt((1.-rstate[0])*(1.-rstate[1])) ; ret[0][2] = 0 ;
	ret[1][0] = m[1][0]*sqrt((1.-rstate[0])*(1.-rstate[1]))  ; ret[1][1] = m[1][1]*(1.-rstate[1])  ; ret[1][2] = 0 ;
	ret[2][0] = 0 ; ret[2][1] = 0 ; ret[2][2] = m[2][2]*sqrt((1.-rstate[0])*(1.-rstate[1])) ;
	return ret ;

}


bool AnisotropicLinearDamage::fractured() const
{
	if (fraction < 0)
		return false ;
	
// 	if(inTension)
// 		if(std::min(getState()[2], std::min(getState()[1],getState()[3])) >= secondaryThresholdDamageDensity)
// 			return true ;
// 		
// 	if(inCompression)
// 		if(getState()[0] >= thresholdDamageDensity)
// 			return true ;
		
// 	std::cout << std::max(tensionDamage, compressionDamage) <<  " " << thresholdDamageDensity/**fraction*/ << std::endl ;
	return getState().min() >= thresholdDamageDensity ;
}

AnisotropicLinearDamage::~AnisotropicLinearDamage()
{
}


}
