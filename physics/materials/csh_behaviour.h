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
    CSHBehaviour(CSHType type, Function ageingFunction = Function("1"), double E=16e9, double nu=0.25,  SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;
    CSHBehaviour(CSHType type, const  std::vector<double> & densities, const std::vector<double> & times, double E=16e9, double nu=0.25,  SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;
    double currentTime ;
    CSHType CshType ;
    Function ageingFunction ;
    bool function ;
    std::vector<double> densities ;
    std::vector<double> times ;

    virtual Form * getCopy() const ;
    virtual void step(double timestep, ElementState& s, double maxScore);

} ;

} 
#endif // PASTE_BEHAVIOUR
