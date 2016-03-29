// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CONCRETE_BEHAVIOUR_H
#define CONCRETE_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"
#include "../fracturecriteria/mcft.h"

namespace Amie
{

/*PARSE ConcreteBehaviour Form 
    @value[young_modulus] // value of the Young modulus
    @value[poisson_ratio] // value of the Poisson ratio
    @value[compressive_strength] // value of the tensile strength of the material
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
    @string<RedistributionType>[redistribution] UPPER_BOUND // redistribution of the rebars
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
 */
struct ConcreteBehaviour : public WeibullDistributedStiffness
{
    std::vector<std::pair<double, double> > rebarLocationsAndDiameters ;
    RedistributionType rtype ;
    ConcreteBehaviour(double E=37e9, double nu=0.3, double compressive = -37e6, planeType pt = PLANE_STRESS, RedistributionType rtype = UPPER_BOUND,  SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;

    virtual Form * getCopy() const ;

} ;

}

#endif // PASTE_BEHAVIOUR
