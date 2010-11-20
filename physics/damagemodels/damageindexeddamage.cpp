//
// C++ Implementation: isotropiclineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "damageindexeddamage.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"

using namespace Mu ;

IndexedLinearDamage::IndexedLinearDamage(int numDof, double dcost, FractureCriterion * e) : DamageModel(e->getMaterialCharacteristicRadius()), dcost(dcost),e(e)
{
	state.resize(1, 0.);
	fixedDamage.resize(1, 0.);
	isNull = false ;
	thresholdDamageDensity = 1 ;
	this->dcost = 200 ;
	e->setEnergyIndexed(true);
}

const Vector & IndexedLinearDamage::damageState() const
{
	return state ;
}

Vector & IndexedLinearDamage::damageState()
{
	return state ;
}


void IndexedLinearDamage::step(ElementState & s)
{
		double volume = 0 ;
		if(!e->getCache().empty())
		{
			volume = s.getParent()->area() ;
		}
		else if(!e->getCache3d().empty())
		{
			volume = s.getParent()->volume() ;
		}
		
		double dd = 1e-2 ;
		if((1.-state[0])-1e-2 < POINT_TOLERANCE)
			dd = -1e-2 ;
		std::pair<double, double> ener_delta = e->getDeltaEnergyDeltaCriterion(s,dd) ;
		
		//ener_delta.first*delta_d+ener_delta.second*delta_d*dcost = currentEnergy-previousEnergy ;
		//delta_d = (ener_delta.first+ener_delta.second*dcost)/(currentEnergy-previousEnergy)
		double delta_d =e->getDeltaEnergyAtState()/(volume*dcost +ener_delta.first) ;
		if(delta_d < -ener_delta.second*state[0])
			delta_d = -ener_delta.second*state[0] ;
		if(state[0]+delta_d >= 1)
			state[0] = 1 ;
		else if(state[0]+delta_d < fixedDamage[0])
			state[0] = fixedDamage[0] ;
		else
			state[0] += delta_d ;
		
		std::cout << ener_delta.second << "  " << delta_d << std::endl ;
		
// 		std::cout << " * " << delta_d << "  " << ener_delta.second*dcost*volume << "  " << ener_delta.first<< std::endl ;
// 		std::cout << " / " << delta_d << "  " << ener_delta.second*dcost/volume << "  " << ener_delta.first<< std::endl ;
// 		std::cout << " . " << delta_d << "  " << e->getDeltaEnergyAtState()<< "  " << ener_delta.first<< std::endl ;
}

void IndexedLinearDamage::artificialDamageStep(double d)
{
	state[0] = std::min(state[0]+d,thresholdDamageDensity/fraction+POINT_TOLERANCE) ;
}


Matrix IndexedLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0. ;
	return ret*(1.-state[0]) ;
}


Matrix IndexedLinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0. ;
	return ret*(1.-previousstate[0]) ;
}

bool IndexedLinearDamage::fractured() const 
{
// 	if(fraction < 0)
		return false ;
// 	return state[0] >= thresholdDamageDensity/fraction ;
}


IndexedLinearDamage::~IndexedLinearDamage()
{
}


