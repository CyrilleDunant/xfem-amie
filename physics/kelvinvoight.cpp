//
// C++ Implementation: kelvinvoight
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "kelvinvoight.h"

using namespace Mu ;

KelvinVoight::KelvinVoight(const Matrix & rig, const Matrix & e) : LinearForm(rig, false, false, rig.numRows()/3+1), eta(e)
{
} ;

KelvinVoight::~KelvinVoight() { } ;

Matrix KelvinVoight::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
	
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) 
    +      VirtualMachine().ieval(Differential(TIME_VARIABLE)*(Gradient(p_i) * eta * Gradient(p_j, true)), e,v);
}

Matrix KelvinVoight::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);

	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v) 
	+      VirtualMachine().ieval(GradientDot(p_i) * eta * GradientDot(p_j, true), gp, Jinv,v);
}

bool KelvinVoight::fractured() const
{
	return false ;
}

Form * KelvinVoight::getCopy() const 
{
	return new KelvinVoight(*this) ;
}

Vector KelvinVoight::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{
	return Vector(0) ;
}

