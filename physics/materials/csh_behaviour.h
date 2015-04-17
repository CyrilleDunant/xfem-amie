// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2015
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CSH_BEHAVIOUR_H
#define CSH_BEHAVIOUR_H

#include "../stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
typedef enum
{
    INNER_CSH,
    OUTER_CSH
} CSHType ;

struct CSHBehaviour : public Stiffness
{
    // per http://www.sciencedirect.com/science/article/pii/S1359645408008720
    //    0.32    0.40    0.48
    //            1.21
    //h   36
    //l    
bool fromPorosity ;
    CSHBehaviour(CSHType type, const  std::vector<double> & densities, const std::vector<double> & times, bool fromPorosity= false, double E=24e9, double nu=0.25, SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;
    double currentTime ;
    CSHType CshType ;
//     Function ageingFunction ;
//     bool function ;
    std::vector<double> densities ;
    std::vector<double> times ;

    virtual Form * getCopy() const ;
    virtual void step(double timestep, ElementState& s, double maxScore);

} ;

} 
#endif // PASTE_BEHAVIOUR
