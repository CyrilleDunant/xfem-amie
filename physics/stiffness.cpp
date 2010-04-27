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
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
} ;

Stiffness::~Stiffness() { } ;

Matrix Stiffness::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}

void Stiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
// std::cout << "a--" << std::endl ;
// Jinv[0].print() ;
// std::cout << "--b" << std::endl ;
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;
// 	ret.print() ;
}

bool Stiffness::fractured() const
{
	return false ;
}

Form * Stiffness::getCopy() const 
{
	return new Stiffness(*this) ;
}

void Stiffness::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
}


PseudoPlastic::PseudoPlastic(const Mu::Matrix& rig, FractureCriterion* crit, DamageModel * damagemodel): LinearForm(rig, false, true, rig.numRows()/3+1), crit(crit), damagemodel(damagemodel), alpha(1), change(true)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
}

bool PseudoPlastic::changed() const
{
	return change ;
} 

PseudoPlastic::~PseudoPlastic() { } ;

Matrix PseudoPlastic::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	return VirtualMachine().ieval(Gradient(p_i) * (param*alpha) * Gradient(p_j, true), e,v) ;
}

void PseudoPlastic::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
// std::cout << "a--" << std::endl ;
// Jinv[0].print() ;
// std::cout << "--b" << std::endl ;
	vm->ieval(Gradient(p_i) * (param*alpha) * Gradient(p_j, true), gp, Jinv,v, ret) ;
// 	ret.print() ;
}

void PseudoPlastic::step(double timestep, ElementState & currentState)
{
	change = false ;
	double prevalpha = alpha ;
	Vector v = currentState.getPrincipalStresses(currentState.getParent()->getCenter()) ;
	Matrix compliance = inverse3x3Matrix(param*alpha) ;
	Vector eps = compliance*v ;
	
	double score = sqrt(eps[0]*eps[0]+eps[1]*eps[1]+eps[2]*eps[2]) ;
	if(score > 0.02*(1.-0.000001))
		alpha = 0.000001 ;
	else if (score > 0)
	{
		alpha = 1.-(score - 0)/(0.02-0) ;
	}
	else
		alpha = 1. ;
	
	change = std::abs(alpha - prevalpha) > POINT_TOLERANCE ;
}

Matrix PseudoPlastic::getTensor(const Point & p) const
{
	return (param*alpha) ;
}

bool PseudoPlastic::fractured() const
{
	return false ;
}

Form * PseudoPlastic::getCopy() const 
{
	return new PseudoPlastic(*this) ;
}

void PseudoPlastic::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
}


