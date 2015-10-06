//
// C++ Implementation: isotropiclineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "damageindexeddamage.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"

using namespace Amie ;

IndexedLinearDamage::IndexedLinearDamage(int numDof, double dcost, FractureCriterion * e) : dcost(dcost),e(e)
{
	state.resize(1, 0.);
	fixedDamage.resize(1, 0.);
	isNull = false ;
	thresholdDamageDensity = 1 ;
	this->dcost = 200 ;
	e->setEnergyIndexed(true);
}

std::pair<Vector, Vector> IndexedLinearDamage::computeDamageIncrement(ElementState & s)
{
		Vector ret(1) ; ret = 0 ;
		
		double detot = 0 ;
		double deavg = 0 ;
		double vtot = 0 ;
		double remnantEnergy = e->getDeltaEnergyAtState() ;
		double originalEnergy = remnantEnergy ;
		std::vector<DelaunayTriangle *> totry ;
		for(size_t i = 0 ; i < e->getCache().size() ; i++)
		{
// 			if(e->getCache()[i]->getBehaviour()->getFractureCriterion()->met(e->getCache()[i]->getState()))
				totry.push_back(static_cast<DelaunayTriangle *>(e->mesh2d->getInTree(e->getCache()[i]))) ;
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
			if(maxdd + totry[mindeddindex]->getBehaviour()->getDamageModel()->getState()[0] > 1)
				maxdd = 1.-totry[mindeddindex]->getBehaviour()->getDamageModel()->getState()[0] ;
			if(maxdd + totry[mindeddindex]->getBehaviour()->getDamageModel()->getState()[0] < 0)
				maxdd = totry[mindeddindex]->getBehaviour()->getDamageModel()->getState()[0] ;
			
			maxdd *=.01 ;
			if(totry[mindeddindex] == s.getParent())
			{

				ret[0] += maxdd ;
				return std::make_pair(state,ret);
			}
			remnantEnergy += maxdd*totry[mindeddindex]->area()*totry[mindeddindex]->getBehaviour()->getFractureCriterion()->getEnergyDamageDifferential() ;
			totry.erase(totry.begin()+mindeddindex) ;
		}
		
		for(size_t i = 0 ; i < e->getCache().size() ; i++)
		{
			DelaunayTriangle * tri = static_cast<DelaunayTriangle *>( e->mesh2d->getInTree(e->getCache()[i]) ) ;
			double dedd =tri->getBehaviour()->getFractureCriterion()->getEnergyDamageDifferential() ;
			double a = tri->area() ;
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
		if(state[0] + ret[0] + delta_d >= 1)
			ret[0] = 1 ;
		else if(state[0] + ret[0] + delta_d < fixedDamage[0])
			ret[0] = fixedDamage[0] ;
		else
			ret[0] += delta_d ;
		
		return std::make_pair(state,ret) ;

}

void IndexedLinearDamage::computeDelta(ElementState & s)
{
	Vector ret(1) ; ret = 0 ;
	
	double detot = 0 ;
	double deavg = 0 ;
	double vtot = 0 ;
	double remnantEnergy = e->getDeltaEnergyAtState() ;
	double originalEnergy = remnantEnergy ;
	std::vector<DelaunayTriangle *> totry ;
	for(size_t i = 0 ; i < e->getCache().size() ; i++)
	{
		totry.push_back(static_cast<DelaunayTriangle *>(e->mesh2d->getInTree(e->getCache()[i]))) ;
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
		if(maxdd + totry[mindeddindex]->getBehaviour()->getDamageModel()->getState()[0] > 1)
			maxdd = 1.-totry[mindeddindex]->getBehaviour()->getDamageModel()->getState()[0] ;
		if(maxdd + totry[mindeddindex]->getBehaviour()->getDamageModel()->getState()[0] < 0)
			maxdd = totry[mindeddindex]->getBehaviour()->getDamageModel()->getState()[0] ;
		
		maxdd *=.01 ;
		if(totry[mindeddindex] == s.getParent())
		{

			ret[0] += maxdd ;
			delta = (ret-state).max() ;
			return ;
		}
		remnantEnergy += maxdd*totry[mindeddindex]->area()*totry[mindeddindex]->getBehaviour()->getFractureCriterion()->getEnergyDamageDifferential() ;
		totry.erase(totry.begin()+mindeddindex) ;
	}
	
	for(size_t i = 0 ; i < e->getCache().size() ; i++)
	{
		DelaunayTriangle * tri = static_cast<DelaunayTriangle *>( e->mesh2d->getInTree(e->getCache()[i]) ) ;
		double dedd =tri->getBehaviour()->getFractureCriterion()->getEnergyDamageDifferential() ;
		double a = tri->area() ;
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
		if(state[0] + ret[0] + delta_d >= 1)
			ret[0] = 1 ;
		else if(state[0] + ret[0] + delta_d < fixedDamage[0])
			ret[0] = fixedDamage[0] ;
		else
			ret[0] += delta_d ;
		
		std::cout << delta_d << std::endl ;
		delta = (ret-state).max() ;
}

Matrix IndexedLinearDamage::apply(const Matrix & m, const Point & p , const IntegrableEntity * e , int g ) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*1e-6 ;
	return ret*(1.-state[0]) ;
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

DamageModel * IndexedLinearDamage::getCopy() const
{
    IndexedLinearDamage * ret = new IndexedLinearDamage(1,dcost, e) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}

