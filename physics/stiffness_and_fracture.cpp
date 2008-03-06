//
// C++ Interface: stiffness_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_and_fracture.h"
#include "../delaunay.h"
#include "ruptureenergy.h"


using namespace Mu ;


StiffnessAndFracture::StiffnessAndFracture(const Matrix & rig, FractureCriterion * crit) : LinearForm(rig, false, true, rig.numRows()/3+1) 
{
	criterion = crit ;
	frac = false ;
	init = param[0][0] ;
	change  = false ;
	previousDamage = 0 ;
	damage = 0 ;
} ;

StiffnessAndFracture::~StiffnessAndFracture() 
{ 
	delete criterion ;
} ;

Matrix StiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36 )
		v.push_back(ZETA);
	
	VirtualMachine vm ;
	return vm.ieval(Gradient(p_i) * (param*(1.-damage)) * Gradient(p_j, true), e,v) ;
}

Matrix StiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	
	return VirtualMachine().ieval(Gradient(p_i) * (param*(1.-damage)) * Gradient(p_j, true), gp, Jinv,v) ;
}


void StiffnessAndFracture::stepBack()
{
	damage = previousDamage ;
}

void StiffnessAndFracture::step(double timestep, ElementState & currentState) 
{
	change = false ;

	if(!frac && criterion->met(currentState) )
	{
		previousDamage = damage ;
		
		damage += .1/**currentState.getParent()->area()*1000000.*/ ;
		change = true ;
		if(damage > .5)
		{
			frac = true ;
			damage = .999 ;
// 			param[0][1] = 0 ;param[0][1] = 0 ;
// 			param[2][2] *= 0.0001 ;
// 			this->type = VOID_BEHAVIOUR ;
		}
	}

}

bool StiffnessAndFracture::changed() const
{
	return change ;
} 

bool StiffnessAndFracture::fractured() const
{
	return frac;
}

Form * StiffnessAndFracture::getCopy() const 
{
	StiffnessAndFracture * copy = new StiffnessAndFracture(param, criterion->getCopy()) ;
	copy->damage = damage ;
	return copy ;
}

Vector StiffnessAndFracture::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{
	return Vector(0) ;
}

Matrix StiffnessAndFracture::getTensor(const Point & p) const
{
	return param*(1.-damage) ;
}

