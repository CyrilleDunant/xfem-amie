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

Material Stiffness::toMaterial()
{
	Material mat(getTensor(Point(0,0))) ;
	return mat ;
}


PseudoPlastic::PseudoPlastic(const Mu::Matrix& rig, double limitStrain, double radius): LinearForm(rig, false, true, rig.numRows()/3+1), limitStrain(limitStrain), radius(radius), alpha(1), change(true)
{

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

PseudoPlastic::~PseudoPlastic() { } ;


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
	if(timestep > POINT_TOLERANCE)
		fixLastDamage() ;
	frac = fixedfrac ;
	change = false ;
	double lastalpha = alpha ;
	if(cache.empty())
	{
		Circle c(2.*radius, currentState.getParent()->getCenter()) ;
		cache = currentState.getParent()->get2DMesh()->getConflictingElements(&c) ;
	}
	double area = 0 ;
	double str = 0 ;
	double fact = 0 ;
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		if( cache[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			double d =  exp(-squareDist2D(currentState.getParent()->getCenter(), cache[i]->getCenter())/(radius*radius)) ;
			double a = cache[i]->area() ;
			str += cache[i]->getState().getVonMisesStrain(cache[i]->getCenter())*a*d ;
			area += a ;
			fact+=d*a ;
		}
	}
	
	double maxStrain = std::abs(str)/fact ;
	
	if(maxStrain > limitStrain)
	{
		currentState.getParent()->behaviourUpdated = true ;
		alpha = std::min(limitStrain/maxStrain, lastDamage) ;
		change = std::abs(lastalpha-alpha) > 1e-4 ;
	}
}

Matrix PseudoPlastic::getTensor(const Point & p) const
{
	return (param*alpha) ;
}

Matrix PseudoPlastic::getPreviousTensor(const Point & p) const
{
	return (param*lastDamage) ;
}

FractureCriterion * PseudoPlastic::getFractureCriterion() const
{
	return NULL ;
}

bool PseudoPlastic::fractured() const
{
	return frac ;
}

Form * PseudoPlastic::getCopy() const 
{
	return new PseudoPlastic(param, limitStrain, radius) ;
	
}


