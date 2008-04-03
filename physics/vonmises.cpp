//
// C++ Implementation: vonmises
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "vonmises.h"
#include "../delaunay.h"
namespace Mu {

VonMises::VonMises(double thresh) : threshold(thresh)
{
}


VonMises::~VonMises()
{
}

double VonMises::grade(const ElementState &s) const
{
	double maxStress = s.getMaximumVonMisesStress() ;

	if(maxStress > threshold )
	{
		return 1. - std::abs(threshold/maxStress) ;
	}
	else 
	{
		return 0 ;
	}
}

bool VonMises::met(const ElementState & s)
{
	DelaunayTriangle * tested = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	if(tested)
	{
		std::vector<DelaunayTriangle *> neighbourhood ;
		for(size_t i = 0 ; i< tested->neighbourhood.size() ; i++)
			neighbourhood.push_back(tested->getNeighbourhood(i)) ;
		
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
	else
	{
		std::cout << " criterion not implemented for this kind of element" << std::endl ;
		return false ;
	}
	
	//shut up the compiler
	return false ;
	
}

FractureCriterion * VonMises::getCopy() const
{
	return new VonMises(*this) ;
}

}
