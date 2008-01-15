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

bool MohrCoulomb::met(const ElementState & s) const
{
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;
	if(testedTri)
	{
		std::set<DelaunayTriangle *> neighbourhood ;
		std::vector<DelaunayTriangle *> neighbours = testedTri->neighbourhood ;
		for(size_t i = 0 ; i < neighbours.size() ; i++)
		{
			for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
			{
				if(neighbours[i]->neighbourhood[j] != testedTri 
				   && !neighbours[i]->neighbourhood[j]->getBehaviour()->fractured())
					neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
			}
		}

		double maxNeighbourhoodStress = 0 ;
		double minNeighbourhoodStress = 0 ;
		if(!neighbourhood.empty())
		{
			for(std::set<DelaunayTriangle *>::const_iterator i = neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
			{
				Vector pstress = (*i)->getState().getPrincipalStresses((*i)->getCenter()) ;
				double maxStress = pstress.max() ;
				double minStress = pstress.min() ;
				if(maxStress > maxNeighbourhoodStress)
					maxNeighbourhoodStress = maxStress ;
				if(minStress < minNeighbourhoodStress)
					minNeighbourhoodStress = minStress ;
			}
		}
		
		Vector pstress = s.getPrincipalStresses(Point(1./3., 1./3.), true) ;
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
	else if(testedHex)
	{
		std::set<HexahedralElement *> neighbourhood ;
		std::vector<HexahedralElement *> neighbours = testedHex->neighbourhood ;
		for(size_t i = 0 ; i < neighbours.size() ; i++)
		{
			for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
			{
				if(neighbours[i]->neighbourhood[j] != testedHex 
				   && !neighbours[i]->neighbourhood[j]->getBehaviour()->fractured())
					neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
			}
		}

		double maxNeighbourhoodStress = 0 ;
		double minNeighbourhoodStress = 0 ;
		if(!neighbourhood.empty())
		{
			for(std::set<HexahedralElement *>::const_iterator i = neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
			{
				Vector pstress = (*i)->getState().getPrincipalStresses((*i)->getCenter()) ;
				double maxStress = pstress.max() ;
				double minStress = pstress.min() ;
				if(maxStress > maxNeighbourhoodStress)
					maxNeighbourhoodStress = maxStress ;
				if(minStress < minNeighbourhoodStress)
					minNeighbourhoodStress = minStress ;
			}
		}
		
		Vector pstress = s.getPrincipalStresses(Point(1./3., 1./3.), true) ;
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

FractureCriterion * MohrCoulomb::getCopy() const
{
	return new MohrCoulomb(*this) ;
}

}
