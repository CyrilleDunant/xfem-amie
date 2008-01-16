//
// C++ Implementation: stiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness.h"

using namespace Mu ;

Stiffness::Stiffness(const Matrix & rig) : LinearForm(rig, false, false, rig.numRows()/3+1) 
{
} ;

Stiffness::~Stiffness() { } ;

Matrix Stiffness::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}

Matrix Stiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);

	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v) ;
}

bool Stiffness::fractured() const
{
	return false ;
}

Form * Stiffness::getCopy() const 
{
	return new Stiffness(*this) ;
}

Vector Stiffness::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{
	return Vector(0) ;
}

