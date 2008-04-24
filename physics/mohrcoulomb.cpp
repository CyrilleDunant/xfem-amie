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

double MohrCoulomb::grade(const ElementState &s) const 
{
	if(s.getParent()->getBehaviour()->fractured())
		return 0 ;

	Vector pstress = s.getPrincipalStresses(Point(1./3., 1./3.), true) ;
	double maxStress = pstress.max();
	double minStress = pstress.min();
	if( maxStress > upVal )
	{
		return 1. - std::abs(upVal/maxStress) ;
	}
		
	if( minStress < downVal )
	{
		return 1. - std::abs(downVal/minStress) ;
	}

	return 0 ;
}

bool MohrCoulomb::met(const ElementState & s) 
{
	if(s.getParent()->getBehaviour()->fractured())
		return false ;

	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;
	if(testedTri)
	{
		double score = grade(s) ;
		
		if (score == 0)
			return false ;
		
		if(cache.empty())
		{
			Circle epsilon(0.01,testedTri->getCenter()) ;
			cache = testedTri->tree->conflicts(&epsilon);
		}

		double maxNeighbourhoodScore = 0 ;
		if(!cache.empty())
		{
			for(size_t i = 0 ; i< cache.size() ; i++)
			{
				if(cache[i]->getBehaviour()->getFractureCriterion())
					maxNeighbourhoodScore = std::max(maxNeighbourhoodScore,
					                                 cache[i]->getBehaviour()->getFractureCriterion()->grade(cache[i]->getState())) ;
			}
		}
		
		
		if(score >= maxNeighbourhoodScore)
		{
			return true ;
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
				   && !neighbours[i]->neighbourhood[j]->getBehaviour()->fractured())
					neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
			}
		}

		double maxNeighbourhoodScore = 0 ;
		if(!neighbourhood.empty())
		{
			for(std::set<HexahedralElement *>::iterator i= neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
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
	
	//shut up the compiler
	return false ;
}

FractureCriterion * MohrCoulomb::getCopy() const
{
	return new MohrCoulomb(*this) ;
}

}
