// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "csh_behaviour.h"
#include "../homogenization/composite.h"

namespace Amie
{

CSHBehaviour::CSHBehaviour(CSHType type, const  std::vector<double> & densities, const std::vector<double> & times, bool fromPorosity, double E, double nu,  SpaceDimensionality dim ) : Stiffness(Tensor::cauchyGreen(std::make_pair(E,nu), true,dim)),fromPorosity(fromPorosity), currentTime(0), CshType(type), densities(densities), times(times)
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
        if(!fromPorosity)
        {
            double fac = 1. ;
            if(!densities.empty())
            {
                size_t index = 0 ;
                fac = densities[index] ;
                while(times[index] < currentTime && index < times.size())
                    index++ ;
                
                if(index >= times.size() || densities.size() == 1)
                    fac = densities.back() ;
                else 
                {
                    double wb = times[index]-currentTime ;
                    double wa = currentTime-times[index-1] ;
                    fac = (densities[index-1]*wb +  densities[index]*wa)/(wa+wb) ;
                }
            }
                
            return new Stiffness(param*fac) ;
        }
        else
        {
            double phi = 0. ;
            if(!densities.empty())
            {
                size_t index = 0 ;
                phi = densities[index] ;
                while(times[index] < currentTime && index < times.size())
                    index++ ;
                
                if(index >= times.size() || densities.size() == 1)
                    phi = densities.back() ;
                else 
                {
                    double wb = times[index]-currentTime ;
                    double wa = currentTime-times[index-1] ;
                    phi = (densities[index-1]*wb +  densities[index]*wa)/(wa+wb) ;
                }
            }
            
            Phase porosity(new Stiffness(Tensor::cauchyGreen(std::make_pair(0,.4997), true,SPACE_THREE_DIMENSIONAL)), phi) ;
            Phase csh(new Stiffness(Tensor::cauchyGreen(std::make_pair(1,.25), true,SPACE_THREE_DIMENSIONAL)), 1.-phi) ;
            BiphasicSelfConsistentComposite sc(porosity,csh) ;
            double fac = sc.getBehaviour()->getTensor(Point())[0][0]/param[0][0]  ;
            return new Stiffness(param*fac) ;
        }
            

    }
    return new DerivedStiffness(param) ;
}


}
