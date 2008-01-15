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

bool VonMises::met(const ElementState & s) const
{
	DelaunayTriangle * tested = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	if(tested)
	{
		std::vector<DelaunayTriangle *> neighbourhood = tested->neighbourhood ;
		
		double maxNeighbourhoodStress = 0 ;
		if(!neighbourhood.empty())
		{
			for(size_t i = 0 ; i< neighbourhood.size() ; i++)
			{
				double maxStress =  neighbourhood[i]->getState().getMaximumVonMisesStress() ;
				if(maxStress > maxNeighbourhoodStress)
					maxNeighbourhoodStress = maxStress ;

			}
		}
		
		double maxStress = s.getMaximumVonMisesStress() ;

		if(maxStress > threshold )
		{
			if(maxNeighbourhoodStress < maxStress)
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

FractureCriterion * VonMises::getCopy() const
{
	return new VonMises(*this) ;
}

}
