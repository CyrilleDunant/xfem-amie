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
	Material mat(param) ;
	return mat ;
}


PseudoPlastic::PseudoPlastic(const Mu::Matrix& rig, FractureCriterion* crit, DamageModel * damagemodel): LinearForm(rig, false, true, rig.numRows()/3+1), crit(crit), damagemodel(damagemodel), alpha(1), change(true)
{
	lastCritUp = dynamic_cast<MohrCoulomb *>(crit)->upVal ;
	lastCritDown = dynamic_cast<MohrCoulomb *>(crit)->downVal ;
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
	frac = fixedfrac ;
	change = false ;
	Vector p(1) ; p[0] = 1.-lastDamage ;
	if(crit->met(currentState))
	{
		damagemodel->damageState() = p ;
		currentState.getParent()->behaviourUpdated = true ;
// 		double talpha = lastDamage ;
// 		double balpha = 0.00001 ;
		while(crit->grade(currentState) > 0)
		{
// 			alpha = (talpha+balpha)*.5 ;
// 			if(crit->grade(currentState) < 0)
// 				balpha = alpha ;
// 			else
// 				talpha = alpha ;
			double deltaState = damagemodel->damageState()[0] ;
			damagemodel->step(currentState);
			
			if(crit->grade(currentState) < -1e-6)
			{
				while(crit->grade(currentState) < -1e-6)
				{
					damagemodel->state[0]+= 5e-7 ;
					alpha = 1.-damagemodel->damageState()[0] ;
				}
			}
			deltaState = damagemodel->damageState()[0]-deltaState ;
			
// 			DelaunayTriangle * self = dynamic_cast<DelaunayTriangle *>(currentState.getParent()) ;
// 			Circle c(damagemodel->getCharacteristicRadius()*2.,  self->getCenter()) ;
// 			std::vector<DelaunayTriangle *> neighbourhood = self->tree->getConflictingElements(&c) ;
// 			for(int i = 0 ; i < neighbourhood.size() ; i++)
// 			{
// 				if(dynamic_cast<PseudoPlastic *>(neighbourhood[i]->getBehaviour()))
// 				{
// 					double d = dist(self->getCenter(), neighbourhood[i]->getCenter()) ;
// 					PseudoPlastic * psp = dynamic_cast<PseudoPlastic *>(neighbourhood[i]->getBehaviour()) ;
// 					psp->damagemodel->state[0] += deltaState*exp(-d*d/(.25*damagemodel->getCharacteristicRadius()*damagemodel->getCharacteristicRadius())) ;
// 					psp->alpha = 1.-psp->damagemodel->damageState()[0] ;
// 				}
// 			}
			
			alpha = 1.-damagemodel->damageState()[0] ;
			frac = damagemodel->fractured() ;
		}
// 		alpha = balpha ;
// 		std::cout << "c" << std::flush ;
// 		if(frac)
// 			alpha = 0.00001 ;
		change = true ;
	}
	
}

Matrix PseudoPlastic::getTensor(const Point & p) const
{
	return (param*alpha) ;
}

FractureCriterion * PseudoPlastic::getFractureCriterion() const
{
	return crit ;
}

bool PseudoPlastic::fractured() const
{
	return frac ;
}

Form * PseudoPlastic::getCopy() const 
{
	return new PseudoPlastic(param, crit->getCopy(), damagemodel->getCopy()) ;
	
}


