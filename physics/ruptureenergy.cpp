//
// C++ Implementation: ruptureenergy
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ruptureenergy.h"
#include "../delaunay.h"

namespace Mu {

RuptureEnergy::RuptureEnergy(double e)
	: energy(e)
{
}


RuptureEnergy::~RuptureEnergy()
{
}

double RuptureEnergy::grade(const ElementState &s) const
{
	Vector pstress = s.getStress(s.getParent()->getGaussPoints().gaussPoints) ;
	Vector pstrain = s.getStrain(s.getParent()->getGaussPoints().gaussPoints) ;
	
	pstrain *= pstress ;
	
	Vector E(pstrain.size()/3) ;
	
	for(size_t j = 0 ; j < E.size() ; ++j)
	{
		E[j] = pstrain[j*3] + pstrain[j*3+1] + pstrain[j*3+2] ;
	}
	double enr = VirtualMachine().ieval(E, s.getParent())/s.getParent()->area() ;
	if( enr > energy )
	{
		return 1. - std::abs(energy/enr) ;
	}
	else 
	{
		return 0 ;
	}
}

bool RuptureEnergy::met(const ElementState & s) 
{
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;
	if(testedTri)
	{
				Circle epsilon(0.0001,testedTri->getCenter()) ;
		std::vector<DelaunayTriangle *> neighbourhood  = testedTri->tree->conflicts(&epsilon);

		double maxNeighbourhoodScore = 0 ;
		if(!neighbourhood.empty())
		{
			for(size_t i = 0 ; i< neighbourhood.size() ; i++)
			{
				if(neighbourhood[i]->getBehaviour()->getFractureCriterion())
					maxNeighbourhoodScore = std::max(maxNeighbourhoodScore, neighbourhood[i]->getBehaviour()->getFractureCriterion()->grade(neighbourhood[i]->getState())) ;
			}
		}
		
		double score = grade(s) ;
		if( score > 0 )
		{
			if(score > maxNeighbourhoodScore)
			{
				return true ;
			}
		}

		return false ;
		
	}
	else if(testedHex)
	{
		std::set<HexahedralElement *> neighbourhood ;
		std::vector<HexahedralElement *> neighbours = testedHex->neighbourhood ;
		for(size_t i = 0 ; i < neighbours.size() ; i++)
		{
			for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
			{
				if(neighbours[i]->neighbourhood[j] != testedHex 
				   && !neighbours[i]->neighbourhood[j]->getBehaviour()->fractured() 
				   && !(neighbours[i]->neighbourhood[j]->getBehaviour()->type == VOID_BEHAVIOUR))
					neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
			}
		}

		double maxNeighbourhoodScore = 0 ;
		if(!neighbourhood.empty())
		{
			for(std::set<HexahedralElement *>::iterator i = neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
			{
				if((*i)->getBehaviour()->getFractureCriterion())
					maxNeighbourhoodScore = std::max(maxNeighbourhoodScore, (*i)->getBehaviour()->getFractureCriterion()->grade((*i)->getState())) ;
			}
		}
		
		double score = grade(s) ;
		if( score > 0 )
		{
			if(score > maxNeighbourhoodScore)
			{
				return true ;
			}
		}

		return false ;
	}
	else
	{
		std::cout << " criterion not implemented for this kind of element" << std::endl ;
		return false ;
	}
	
	return false ;
}

FractureCriterion * RuptureEnergy::getCopy() const
{
	return new RuptureEnergy(*this) ;
}

}
