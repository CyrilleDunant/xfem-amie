// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "csh_behaviour.h"
#include "../homogenization/composite.h"
#include "../stiffness_with_imposed_stress.h"

namespace Amie
{

CSHBehaviour::CSHBehaviour(CSHType type, const  std::vector<double> & densities, const std::vector<double> & shrinkageStresses, const std::vector<double> & times, double E, double nu,  SpaceDimensionality dim ) : Stiffness(E, nu,dim, PLANE_STRESS, YOUNG_POISSON), currentTime(0), CshType(type), densities(densities), shrinkageStresses(shrinkageStresses), times(times)
{
}

void CSHBehaviour::step(double timestep, ElementState& s, double maxScore)
{
    currentTime += timestep ;
}

Form * CSHBehaviour::getCopy() const 
{
    if(CshType == OUTER_CSH)
    {

	double fac = 1. ;
	double dimensions = (param.size() == 9) ? 2:3 ;
	Vector shrink(0., (size_t)dimensions) ;
	if(!densities.empty())
	{
	    size_t index = 0 ;
	    fac = densities[index] ;
	    while(times[index] < currentTime && index < times.size())
		index++ ;
	    
	    if(index >= times.size() || densities.size() == 1)
	    {
		fac = densities.back() ;
		if(!shrinkageStresses.empty())
		  shrink = shrinkageStresses.back()/dimensions ;
		
	    }
	    else 
	    {
		double wb = times[index]-currentTime ;
		double wa = currentTime-times[index-1] ;
		fac = (densities[index-1]*wb +  densities[index]*wa)/(wa+wb) ;
		if(!shrinkageStresses.empty())
		  shrink = (shrinkageStresses[index-1]*wb +  shrinkageStresses[index]*wa)/(wa+wb) ;
	    }
	}
	
	if(shrinkageStresses.empty())
		return new Stiffness(param*fac) ;	    
	return new StiffnessWithImposedStress(param*fac, shrink) ;

            

    }
    return new DerivedStiffness(param) ;
}


}
