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

bool MaximumStrain::met(const ElementState * s) const
{
	DelaunayTriangle * tested = dynamic_cast<DelaunayTriangle *>(s->getParent()) ;
	if(tested)
	{
		std::vector<DelaunayTriangle *> neighbourhood = tested->neighbourhood ;
		
		double maxNeighbourhoodStrain = 0 ;
		if(!neighbourhood.empty())
		{
			for(size_t i = 0 ; i< neighbourhood.size() ; i++)
			{
				Vector pstrain = neighbourhood[i]->getState().getStrain(neighbourhood[i]->getBoundingPoints()) ;
				double maxStrain = pstrain.max() ;
				if(maxStrain > maxNeighbourhoodStrain)
					maxNeighbourhoodStrain = maxStrain ;
			}
		}
		
		Vector pstrain = s->getStrain(tested->getBoundingPoints()) ;
		double maxStrain = pstrain.max();
		if( maxStrain > upVal )
		{
			if(maxNeighbourhoodStrain < maxStrain)
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
