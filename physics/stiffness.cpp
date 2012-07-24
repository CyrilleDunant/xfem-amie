//
// C++ Implementation: stiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/vonmises.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "fracturecriteria/nonlocalvonmises.h"
#include <valarray>


using namespace Mu ;

Stiffness::Stiffness(const Matrix & rig) : LinearForm(rig, false, false, rig.numRows()/3+1) 
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
} ;

Stiffness::~Stiffness() { } ;

void Stiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
// std::cout << "a--" << std::endl ;
// Jinv[0].print() ;
// std::cout << "--" << std::endl ;
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;
//	ret.print() ;
}

bool Stiffness::fractured() const
{
	return false ;
}

Form * Stiffness::getCopy() const 
{
	return new Stiffness(param) ;
}


PseudoPlastic::PseudoPlastic(const Mu::Matrix& rig, double E, double limitStrain, double radius): LinearForm(rig, false, true, rig.numRows()/3+1), limitStrain(limitStrain), radius(radius), alpha(0), change(true)
{
	stiffness = E ;
	vm = new NonLocalVonMises(limitStrain, E, radius) ;
	vm->setMaterialCharacteristicRadius(radius);
	initialised = false ;
	lastDamage = alpha ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	
	frac = false ;
	fixedfrac = false ;
}

bool PseudoPlastic::changed() const
{
	return change ;
} 

void PseudoPlastic::fixLastDamage()
{
	lastDamage = alpha ;
	fixedfrac = frac ;
}

PseudoPlastic::~PseudoPlastic() { delete vm ;} ;


void PseudoPlastic::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
// std::cout << "a--" << std::endl ;
// Jinv[0].print() ;
// std::cout << "--b" << std::endl ;
	vm->ieval(Gradient(p_i) * (param*(1.-alpha)) * Gradient(p_j, true), gp, Jinv,v, ret) ;
// 	ret.print() ;
}

void PseudoPlastic::step(double timestep, ElementState & currentState, double maxscore)
{
	if(timestep > POINT_TOLERANCE_2D)
	{
		fixLastDamage() ;
	}
	
	frac = fixedfrac ;
	change = false ;
	double lastalpha = alpha ;

	Vector str = vm->smoothedPrincipalStress(currentState) ;
	double maxStress = sqrt( ( ( str[0] - str[1] ) * ( str[0] - str[1] ) + str[0] * str[0] + str[1] * str[1] ) / 2. ) ;
	
	if(maxStress > POINT_TOLERANCE_2D)
	{
		alpha = std::max(1.-(vm->threshold/maxStress)*(1.-alpha), lastDamage) ;
		change = std::abs(alpha-lastalpha) > 1e-6 ;
		currentState.getParent()->behaviourUpdated = change ;
// 		if(change)
// 			std::cout << vm->getScoreAtState() << std::endl ;
	}
}

Matrix PseudoPlastic::getTensor(const Point & p, IntegrableEntity *) const
{
	return (param*(1.-alpha)) ;
}


Matrix PseudoPlastic::getPreviousTensor(const Point & p) const
{
	return (param*lastDamage) ;
}

FractureCriterion * PseudoPlastic::getFractureCriterion() const
{
	return vm ;
}

bool PseudoPlastic::fractured() const
{
	return frac ;
}

Form * PseudoPlastic::getCopy() const 
{
	return new PseudoPlastic(param, stiffness, limitStrain, radius) ;
	
}


