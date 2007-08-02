//
// C++ Implementation: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "mohrcoulomb.h"
#include "../delaunay.h"

namespace Mu {

MohrCoulomb::MohrCoulomb(double up, double down)
	: upVal(up), downVal(down)
{
}


MohrCoulomb::~MohrCoulomb()
{
}

bool MohrCoulomb::met(const ElementState * s) const
{
	DelaunayTriangle * tested = dynamic_cast<DelaunayTriangle *>(s->getParent()) ;
	if(tested)
	{
		std::vector<DelaunayTriangle *> neighbourhood ;
		std::vector<DelaunayTriangle *> neighbours = tested->neighbourhood ;
		for(size_t i = 0 ; i < neighbours.size() ; i++)
		{
			for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
			{
				if(std::find(neighbourhood.begin(), 
				             neighbourhood.end(), 
				             neighbours[i]->neighbourhood[j]) == neighbourhood.end() 
				   && neighbours[i]->neighbourhood[j] != tested)
					neighbourhood.push_back(neighbours[i]->neighbourhood[j]) ;
			}
		}

		double maxNeighbourhoodStress = 0 ;
		double minNeighbourhoodStress = 0 ;
		if(!neighbourhood.empty())
		{
			for(size_t i = 0 ; i< neighbourhood.size() ; i++)
			{
				Vector pstress = neighbourhood[i]->getState()->getPrincipalStresses(neighbourhood[i]->getCenter()) ;
				double maxStress = pstress.max() ;
				double minStress = pstress.min() ;
				if(maxStress > maxNeighbourhoodStress)
					maxNeighbourhoodStress = maxStress ;
				if(minStress < minNeighbourhoodStress)
					minNeighbourhoodStress = minStress ;
			}
		}
		
		Vector pstress = s->getPrincipalStresses(tested->getCenter()) ;
		double maxStress = pstress.max();
		double minStress = pstress.min();
		if( maxStress > upVal )
		{
			if(maxNeighbourhoodStress < maxStress)
			{
				return true ;
			}
		}
		
		if( minStress < downVal )
		{
			if(minNeighbourhoodStress > minStress)
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
	
	//shut up the compiler
	return false ;
}

}
