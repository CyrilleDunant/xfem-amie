// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "csh_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../../utilities/random.h"

namespace Amie
{

CSHBehaviour::CSHBehaviour(CSHType type,Function ageingFunction, double E, double nu, SpaceDimensionality dim) : Stiffness(Material::cauchyGreen(std::make_pair(E,nu), true,dim)), currentTime(0), CshType(type), ageingFunction(ageingFunction), function(true)
{
}

CSHBehaviour::CSHBehaviour(CSHType type, const  std::vector<double> & densities, const std::vector<double> & times, double E, double nu,  SpaceDimensionality dim ) : Stiffness(Material::cauchyGreen(std::make_pair(E,nu), true,dim)), currentTime(0), CshType(type), ageingFunction(Function("1")), function(false), densities(densities), times(times)
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
        if(function)
        {
            double fac = std::min(VirtualMachine().eval(ageingFunction, 0., 0., 0., currentTime), 1.) ;
            return new Stiffness(param*fac) ;
        }
        else
        {
            size_t index = 0 ;
            while(times[index] < currentTime && index < times.size())
                index++ ;
            
            double fac = densities[index] ;
            return new Stiffness(param*fac) ;
            
        }
    }
    return new DerivedStiffness(param) ;
}


}
