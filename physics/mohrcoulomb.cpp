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
		Circle epsilon(0.0005,testedTri->getCenter()) ;
		std::vector<DelaunayTriangle *> neighbourhood ;
		std::vector<DelaunayTriangle *> toTry ;
		for(size_t j = 0 ; j < testedTri->neighbourhood.size() ;j++)
			toTry.push_back(testedTri->getNeighbourhood(j)) ;
		std::vector<DelaunayTriangle *> cleanup ;
		int count = toTry.size() ;
		cleanup.push_back(testedTri) ;
		testedTri->visited = true ;
		while(!toTry.empty())
		{
			std::vector<DelaunayTriangle *> newTrianglesToTry ;
			for(size_t i = 0 ; i < toTry.size() ; i++)
			{
				cleanup.push_back(toTry[i]) ;
				toTry[i]->visited = true ;
				if(epsilon.in(toTry[i]->getCenter()) || count > 0)
				{
					count-- ;
					neighbourhood.push_back(toTry[i]) ;
					for(size_t j = 0 ; j < toTry[i]->neighbourhood.size() ;j++)
					{
						if(!toTry[i]->getNeighbourhood(j)->visited)
							newTrianglesToTry.push_back(toTry[i]->getNeighbourhood(j));
					}
				}
			}
			
			std::sort(newTrianglesToTry.begin(), newTrianglesToTry.end()) ;
			std::vector<DelaunayTriangle *>::iterator e = std::unique(newTrianglesToTry.begin(), newTrianglesToTry.end()) ;
			newTrianglesToTry.erase(e, newTrianglesToTry.end()) ;
			toTry = newTrianglesToTry ;
		}
// 		for(size_t i = 0 ; i < neighbours.size() ; i++)
// 		{
// 			for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
// 			{
// 				if(neighbours[i]->getNeighbourhood(j) != testedTri 
// 				   && !neighbours[i]->getNeighbourhood(j)->getBehaviour()
// 				   && !neighbours[i]->getNeighbourhood(j)->getBehaviour()->fractured())
// 					neighbourhood.insert(neighbours[i]->getNeighbourhood(j)) ;
// 			}
// 		}
		
		for(size_t i = 0 ; i < cleanup.size() ;i++)
			cleanup[i]->visited = false ;


		double maxNeighbourhoodStress = 0 ;
		double minNeighbourhoodStress = 0 ;
		if(!neighbourhood.empty())
		{
			for(std::vector<DelaunayTriangle *>::const_iterator i = neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
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
