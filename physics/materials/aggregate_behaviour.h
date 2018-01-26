// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef AGGREGATE_BEHAVIOUR_H
#define AGGREGATE_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{

/*PARSE AggregateBehaviour Form 
    @string<bool>[elastic] FALSE // desactivates the damage if TRUE
    @string<bool>[space_time] FALSE // uses the behaviour in space-time if TRUE
    @value[young_modulus] 59e9 // value of the Young modulus
    @value[poisson_ratio] 0.3 // value of the Poisson ratio
    @value[tensile_strength] 10e6 // value of the tensile strength of the material
    @value[material_characteristic_radius] 0.00025 // radius of the non-local damage model
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
    @value[variability] 0.2 // variability of the mechanical properties
    @value[blocks] 0 // additional ghost blocks for viscoelastic simulations 
 */
struct AggregateBehaviour : public WeibullDistributedStiffness
{
    int freeblocks ;
    bool spaceTime ;
    bool elastic ;

    AggregateBehaviour(bool elastic = false, bool st = false, double E=59e9, double nu=0.3, double up = 30e6, double r = 0.0004, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS,  double var = 0.2, int blocks = 0) ;

    virtual Form * getCopy() const ;

} ;

}

#endif // AGGREGATE_BEHAVIOUR
