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
	previousstate.resize(1, 0.);
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
	previousstate = state ;
		double volume = 0 ;
		if(!e->getCache().empty())
		{
			volume = s.getParent()->area() ;
		}
		else if(!e->getCache3d().empty())
		{
			volume = s.getParent()->volume() ;
		}
		
		double detot = 0 ;
		double deavg = 0 ;
		double vtot = 0 ;
		double remnantEnergy = e->getDeltaEnergyAtState() ;
		double originalEnergy = remnantEnergy ;
// 		std::cout << "originalEnergy : " << originalEnergy << "  "<< std::flush ;
		std::vector<DelaunayTriangle *> totry ;
		for(size_t i = 0 ; i < e->getCache().size() ; i++)
		{
// 			if(e->getCache()[i]->getBehaviour()->getFractureCriterion()->met(e->getCache()[i]->getState()))
				totry.push_back(e->getCache()[i]) ;
		}
		
		while(!totry.empty() && remnantEnergy > 0)
		{
			int mindeddindex = 0 ;
			double mindedd = totry[0]->getBehaviour()->getFractureCriterion()->getEnergyDamageDifferential() ;
			for(size_t i = 0 ; i < totry.size() ; i++)
			{
				double dedd = totry[i]->getBehaviour()->getFractureCriterion()->getEnergyDamageDifferential() ;
				if(dedd < mindedd)
				{
					mindedd = dedd ;
					mindeddindex = i ;
				}
			}
			double dcdd = totry[mindeddindex]->getBehaviour()->getFractureCriterion()->getCriterionDamageDifferential() ;
			double maxdd = std::max(
				-1./(1.-totry[mindeddindex]->getBehaviour()->getFractureCriterion()->getScoreAtState())/dcdd,
				-originalEnergy/totry[mindeddindex]->getBehaviour()->getFractureCriterion()->getEnergyDamageDifferential() 
 														) ;
			if(std::abs(dcdd) < POINT_TOLERANCE)
				maxdd = 0 ;
			if(maxdd + totry[mindeddindex]->getBehaviour()->getDamageModel()->state[0] > 1)
				maxdd = 1.-totry[mindeddindex]->getBehaviour()->getDamageModel()->state[0] ;
			if(maxdd + totry[mindeddindex]->getBehaviour()->getDamageModel()->state[0] < 0)
				maxdd = totry[mindeddindex]->getBehaviour()->getDamageModel()->state[0] ;
			
			maxdd *=.01 ;
			if(totry[mindeddindex] == s.getParent())
			{
// 				if(std::abs(maxdd) > 1e-7)
// 					std::cout << maxdd << std::endl ;
				state[0] += maxdd ;
				return ;
			}
			remnantEnergy += maxdd*totry[mindeddindex]->area()*totry[mindeddindex]->getBehaviour()->getFractureCriterion()->getEnergyDamageDifferential() ;
			totry.erase(totry.begin()+mindeddindex) ;
		}
		
		for(size_t i = 0 ; i < e->getCache().size() ; i++)
		{
			double dedd = e->getCache()[i]->getBehaviour()->getFractureCriterion()->getEnergyDamageDifferential() ;
			double a = e->getCache()[i]->area() ;
			deavg += dedd*a ;
			vtot += a ;
			detot += dedd*dedd ;
		}
		detot = sqrt(detot) ;
		deavg /= vtot ;
		if(detot < POINT_TOLERANCE)
			detot = 1 ;

// 		std::pair<double, double> ener_delta = e->getDeltaEnergyDeltaCriterion(s,dd) ;
		
		//ener_delta.first*delta_d+ener_delta.second*delta_d*dcost = currentEnergy-previousEnergy ;
		//delta_d = (ener_delta.first+ener_delta.second*dcost)/(currentEnergy-previousEnergy)
		double delta_d = -e->getEnergyDamageDifferential()/detot ; /*/(dcost*(e->getScoreAtState()+ener_delta.second)*.5 +ener_delta.first) */;
// 		delta_d *= volume ;
// 		if(delta_d < -ener_delta.second+state[0])
// 			delta_d = -ener_delta.second+state[0] ;
		if(state[0]+delta_d >= 1)
			state[0] = 1 ;
		else if(state[0]+delta_d < fixedDamage[0])
			state[0] = fixedDamage[0] ;
		else
			state[0] += delta_d ;
		
		std::cout << delta_d << std::endl ;
		
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


