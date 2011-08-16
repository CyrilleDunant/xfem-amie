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

PlasticStrain::PlasticStrain() : previousImposedStrain(0.,3), imposedStrain(0.,3)
{
	getState(true).resize(1, 0.);
	getPreviousState().resize(1, 0.);
	isNull = false ;
	es = NULL ;
	v.push_back(XI);
	v.push_back(ETA);
	param = NULL ;
}

std::pair<Vector, Vector> PlasticStrain::computeDamageIncrement(ElementState & s)
{
	Vector ret = s.getStrain(s.getParent()->getCenter());
	if(ret.size() > 3 && v.size() == 2)
		v.push_back(ZETA);
	if(!param)
		param = new Matrix(s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;
	if(!es)
		es =&s ;
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() && s.getParent()->getBehaviour()->getFractureCriterion()->met())
	{
		if(v.size() == 2 && !previousImposedStrain.size())
			previousImposedStrain.resize(3, 0.) ;
		
		if(v.size() == 3 && !previousImposedStrain.size())
			previousImposedStrain.resize(6, 0.) ;
		
		if(v.size() == 2 && !imposedStrain.size())
			imposedStrain.resize(3, 0.) ;
		
		if(v.size() == 3 && !imposedStrain.size())
			imposedStrain.resize(6, 0.) ;
		
		imposedStrain = -(s.getStress(s.getParent()->getCenter())-lastStress)*20. ;
		//we only use the deviatoric stress
		double diag = (imposedStrain[0]+imposedStrain[1])*.5 ;
		if(v.size() == 3)
			diag = (imposedStrain[0]+imposedStrain[1]+imposedStrain[2])*.333333333333333333333333333 ;
		imposedStrain[0] -= diag ;
		imposedStrain[1] -= diag ;
		if(v.size() == 3)
			imposedStrain[2] -= diag ;
		
	}

	return std::make_pair( Vector(0., 1), Vector(1., 1)) ;

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
	double m = 0 ;
	
	std::vector<BoundaryCondition * > ret ;
	if(!param || !imposedStrain.size())
		return ret ;
	Vector f = VirtualMachine().ieval(Gradient(p_i) * (imposedStrain*state[0]+previousImposedStrain), gp, Jinv,v) ;
// 	if(state[0] > POINT_TOLERANCE_2D)
// 	{
// 		std::cout << "strain = "<< (imposedStrain*state[0]+previousImposedStrain)[0] << "  " << (imposedStrain*state[0]+previousImposedStrain)[1] << std::endl ;
// 		std::cout << "delta = "<< f[0] << "  " << f[1] << std::endl ;
// 	}

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
	if(v.size() == 2 && !param)
		return Vector(0., 3) ;
	if(v.size() == 3 && !param)
		return Vector(0., 6) ;
		
	return (imposedStrain*state[0]+previousImposedStrain) ;
}

bool PlasticStrain::fractured() const 
{
	if(fraction < 0)
		return false ;
	return state.max() >= 1. ;
}

void PlasticStrain::postProcess()
{
	
	if(converged)
	{
		if(v.size() == 2 && lastStress.size() == 0)
			lastStress.resize(3, 0.);
		if(v.size() == 3 && lastStress.size() == 0)
			lastStress.resize(6, 0.);
		
		lastStress = es->getStress(es->getParent()->getCenter()) ;
	}
	
	if(converged && state[0] > POINT_TOLERANCE_2D)
	{

		previousImposedStrain += imposedStrain*state[0] ;
		
		state[0] = 0;
		imposedStrain = 0 ;
		 

		wasBroken = false ;
	}

}

PlasticStrain::~PlasticStrain()
{
	delete param ;
}


}
