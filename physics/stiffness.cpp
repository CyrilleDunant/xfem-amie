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
	lastDamage = alpha ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	
	frac = false ;
}

bool PseudoPlastic::changed() const
{
	return change ;
} 

void PseudoPlastic::fixLastDamage()
{
	lastDamage = alpha ;
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

	double volume ;
	if(currentState.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
		volume = currentState.getParent()->area() ;
	else
		volume = currentState.getParent()->volume() ;
		
	double charVolume ;
	if(currentState.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
		charVolume = M_PI*damagemodel->getCharacteristicRadius()*damagemodel->getCharacteristicRadius() ;
	else
		charVolume = 4./3.*M_PI*damagemodel->getCharacteristicRadius()*damagemodel->getCharacteristicRadius()*damagemodel->getCharacteristicRadius() ;
	double fraction = volume/charVolume ;

	
	double score = crit->grade(currentState) ;
	double palpha = lastDamage ;
	double prevalpha = alpha ;
	double dalpha = 1.-damagemodel->getThresholdDamageDensity()/fraction ;
	dalpha = std::max(dalpha, 0.00001) ;
// 	std::cout << dalpha << "  "<< alpha << "  "<< palpha << "  "<< score << std::endl ;
	if(std::abs(palpha-dalpha) < POINT_TOLERANCE)
	{
		alpha = 0.00001 ;
	}
	else if( score > 1e-5)
	{
		while( std::abs(score) > 1e-5)
		{
			alpha = (palpha+dalpha)*.5 ;
			score = crit->grade(currentState) ;
			
			if(score < 0)
			{
				dalpha = alpha ;
			}
			else
			{
				palpha = alpha ;
			}
			
		}
	}
	
	if(alpha < 1.-damagemodel->getThresholdDamageDensity()/fraction)
	{
		alpha = 0.00001 ;
		frac = true ;
	}
	change = std::abs(alpha - prevalpha) > POINT_TOLERANCE ;
}

Matrix PseudoPlastic::getTensor(const Point & p) const
{
	return (param*alpha) ;
}

bool PseudoPlastic::fractured() const
{
	return frac ;
}

Form * PseudoPlastic::getCopy() const 
{
	return new PseudoPlastic(*this) ;
}

void PseudoPlastic::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
}


