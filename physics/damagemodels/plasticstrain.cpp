//
// C++ Implementation: isotropiclineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "plasticstrain.h"
#include "../../features/boundarycondition.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu {

PlasticStrain::PlasticStrain() 
{
	getState(true).resize(2, 0.);
	getPreviousState().resize(2, 0.);
	isNull = false ;
	v.push_back(XI);
	v.push_back(ETA);
	param = NULL ;
}

Vector PlasticStrain::computeDamageIncrement(ElementState & s)
{
	Vector ret(0., 2) ;
	ret[0] = 1. ;
	if(!param)
		param = new Matrix(s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;
	if(s.getDeltaTime() > POINT_TOLERANCE_2D || s.getParent()->getBehaviour()->getFractureCriterion() && s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint())
	{
// 		std::cout << "plouf" << std::endl ;
		Vector strain = s.getStrain(s.getParent()->getCenter()) ;
		strain[2] = 0 ;
		if(imposedStrain.size() != strain.size())
		{
			
			previousImposedStrain.resize(strain.size(), 0.);
			imposedStrain.resize( strain.size(), 0.);
			if(strain.size() > 3)
				v.push_back(ZETA);
		}
		Vector is = state[0]*imposedStrain + previousImposedStrain;
		imposedStrain = (strain-is)*(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()) ;
		previousImposedStrain = is ; //previousImposedStrain+is ;
		
// 		std::cout << previousImposedStrain.max() << "  " << imposedStrain.max() << "   "<< s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() << std::endl ;
		state[1] += state[0] ;
		state[0] = 0 ;
	}
	
	if(smoothState.size() != state.size())
		smoothState.resize(state.size());
	
	smoothState = smoothedState(s) ;
	return ret ;

}

void PlasticStrain::artificialDamageStep(double d)
{
	getState(true)[0] = std::min(getState()[0]+d,thresholdDamageDensity/fraction+POINT_TOLERANCE_2D) ;
}


Matrix PlasticStrain::apply(const Matrix & m) const
{
	return m ;
}


Matrix PlasticStrain::applyPrevious(const Matrix & m) const
{
	return m ;
}

std::vector<BoundaryCondition * > PlasticStrain::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	

	std::vector<BoundaryCondition * > ret ;
	if(!imposedStrain.size() || !param)
		return ret ;
	
	Vector f = VirtualMachine().ieval(Gradient(p_i) * (*param * (previousImposedStrain + (/*2.*smoothState[0]-*/state[0])*imposedStrain)), gp, Jinv,v) ;
	
	if(f.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[1]));
	}
	if(f.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[2]));
	}
	return ret ;
}

Vector PlasticStrain::getImposedStress(const Point & p) const
{
	if(!param)
	{
		if(v.size() == 2)
			return Vector (0., 3) ;
		return Vector (0., 6) ;
	}
	return *param*(previousImposedStrain + (/*2.*smoothState[0]-*/state[0])*imposedStrain) ;
}

bool PlasticStrain::fractured() const 
{
	if(fraction < 0)
		return false ;
	return state[0] >= thresholdDamageDensity ;
}


PlasticStrain::~PlasticStrain()
{
	delete param ;
}


}
