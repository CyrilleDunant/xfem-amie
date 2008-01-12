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

bool RuptureEnergy::met(const ElementState & s) const
{
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;
	if(testedTri)
	{
		Circle epsilon(.01, s.getParent()->getCenter()) ;
		std::set<DelaunayTriangle *> neighbourhood ;
		std::vector<DelaunayTriangle *> neighbours = testedTri->neighbourhood ;
		std::vector<DelaunayTriangle *> cleanup ;
		while(!neighbours.empty())
		{
			std::vector<DelaunayTriangle *> totest ;
			for(size_t i = 0 ; i < neighbours.size() ; i++)
			{
				if( !neighbours[i]->visited 
				    && neighbours[i] != testedTri 
					&& (epsilon.in(*neighbours[i]->first) 
				       || epsilon.in(*neighbours[i]->second)
				       || epsilon.in(*neighbours[i]->third) 
				       || neighbours[i]->in(epsilon.getCenter()))
				  )
				{
					cleanup.push_back(neighbours[i]) ;
					neighbours[i]->visited = true ;
					
					if(!neighbours[i]->getBehaviour()->fractured() && !(neighbours[i]->getBehaviour()->type == VOID_BEHAVIOUR))
						neighbourhood.insert(neighbours[i]) ;
					
					totest.insert(totest.end(), neighbours[i]->neighbourhood.begin(), neighbours[i]->neighbourhood.end()) ;
				}
			}
			
			std::sort(totest.begin(), totest.end()) ;
			std::vector<DelaunayTriangle *>::iterator e = std::unique(totest.begin(), totest.end()) ;
			totest.erase(e, totest.end()) ;
			neighbours = totest ;
		}
		
		for(size_t i = 0 ; i < cleanup.size() ; i++)
			cleanup[i]->clearVisited() ;
		
// 		for(size_t i = 0 ; i < neighbours.size() ; i++)
// 		{
// 			for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
// 			{
// 				if(neighbours[i]->neighbourhood[j] != testedTri 
// 				   && !neighbours[i]->neighbourhood[j]->getBehaviour()->fractured() 
// 				   && !(neighbours[i]->neighbourhood[j]->getBehaviour()->type == VOID_BEHAVIOUR))
// 					neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
// 			}
// 		}

		double maxNeighbourhoodEnergy = 0 ;

		if(!neighbourhood.empty())
		{
			for(std::set<DelaunayTriangle *>::const_iterator i = neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
			{
				Vector pstress = (*i)->getState().getStress((*i)->getGaussPoints().gaussPoints) ;
				Vector pstrain = (*i)->getState().getStrain((*i)->getGaussPoints().gaussPoints) ;
				Vector E(pstrain.size()/3) ;
				
				for(size_t j = 0 ; j < E.size() ; ++j)
				{
					E[j] = pstrain[j*3] + pstrain[j*3+1] + pstrain[j*3+2] ;
				}
				
				double enr = VirtualMachine().ieval(E, (*i))/(*i)->area() ;
				
				if(enr > maxNeighbourhoodEnergy)
					maxNeighbourhoodEnergy = enr ;
			}
		}

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
			if(maxNeighbourhoodEnergy < enr)
			{
				return true ;
			}
		}
		
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

		double maxNeighbourhoodEnergy = 0 ;
		if(!neighbourhood.empty())
		{
			for(std::set<HexahedralElement *>::const_iterator i = neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
			{
				Vector pstress = (*i)->getState().getStress((*i)->getGaussPoints().gaussPoints) ;
				Vector pstrain = (*i)->getState().getStrain((*i)->getGaussPoints().gaussPoints) ;
				pstrain *= pstress ;
				Vector E(pstrain.size()/3) ;
				
				
				for(size_t j = 0 ; j < E.size() ; ++j)
				{
					E[j] = pstrain[j*3] + pstrain[j*3+1] + pstrain[j*3+2] ;
				}
				
				double enr = VirtualMachine().ieval(E, (*i))/(*i)->volume() ;
				if(enr > maxNeighbourhoodEnergy)
					maxNeighbourhoodEnergy = enr ;
			}
		}
		Vector pstress = s.getStress(s.getParent()->getGaussPoints().gaussPoints) ;
		Vector pstrain = s.getStrain(s.getParent()->getGaussPoints().gaussPoints) ;
		pstrain *= pstress ;
		Vector E(pstrain.size()/3) ;
		
		for(size_t j = 0 ; j < E.size() ; ++j)
		{
			E[j] = pstrain[j*3] + pstrain[j*3+1] + pstrain[j*3+2] ;
		}
		double enr = VirtualMachine().ieval(E, s.getParent())/s.getParent()->volume() ;
		if( enr > energy )
		{
			if(maxNeighbourhoodEnergy < enr)
			{
				return true ;
			}
		}
	}
	else
	{
		std::cout << " criterion not implemented for this kind of element" << std::endl ;
		return false ;
	}
	
	return false ;
}

}
