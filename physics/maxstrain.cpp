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
#include "maxstrain.h"
#include "../delaunay.h"

namespace Mu {

MaximumStrain::MaximumStrain(double up)
	: upVal(up)
{
}


MaximumStrain::~MaximumStrain()
{
}

double MaximumStrain::grade(const ElementState &s) const
{
	Vector pstrain = s.getStrain(s.getParent()->getBoundingPoints()) ;
	double maxStrain = pstrain.max();
	if(maxStrain > upVal)
		return 1.-std::abs(upVal/maxStrain) ;
	else
		return 0 ;
	
}

bool MaximumStrain::met(const ElementState & s) 
{
	DelaunayTriangle * tested = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	if(tested)
	{
		Circle epsilon(0.0001,tested->getCenter()) ;
		std::vector<DelaunayTriangle *> neighbourhood  = tested->tree->conflicts(&epsilon);

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

FractureCriterion * MaximumStrain::getCopy() const
{
	return new MaximumStrain(*this) ;
}

}
