// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef GEL_BEHAVIOUR_H
#define GEL_BEHAVIOUR_H

#include "../stiffness_with_imposed_deformation.h"
#include "../stiffness_with_imposed_stress.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{

/*PARSE GelBehaviour Form 
    @string<bool>[space_time] FALSE // uses the behaviour in space-time if TRUE
    @value[young_modulus] // value of the Young modulus
    @value[poisson_ratio] // value of the Poisson ratio
    @value[imposed_deformation] // value of the imposed deformation of the material
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
    @value[blocks] 0 // additional ghost blocks for viscoelastic simulations 
 */
struct GelBehaviour : public StiffnessWithImposedStress
{
    bool spaceTime ;
    int freeblocks ;
    
    // variation in density 2.06 to 2.22 of amorphous silica suggests  0.025257
    // volumic expansion is 22%, so linear is 7%,           cube root  0.0685
    // bulteel suggests 0.5 expansion,                             so  0.1447   
    // gel behaviour should be 0.18  
    // Garcia-Diaz and Bulteel "Mechanism of damage for the alkali–silica reaction"  eq. 14 imply 50% volume expansion
    // Leeman values are 8-12 (0.18) GPa "E-modulus of the alkali–silica-reaction product determined by micro-indentation"
    // Moon et al. values are 25 (0.35) GPa "Determination of the elastic properties of amorphous materials: Case study of alkali–silica reaction gel"
    // Chuanlin Hu, Bishnu P. Gautam, Daman K. Panesar, values are 25.7 (+/-7.2) "Nano-mechanical properties of alkali-silica reaction (ASR) productsin concrete measured by nano-indentation"

    GelBehaviour(bool st=false, double E=20e9, double nu=0.35, double alpha= 0.069, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, int blocks = 0) ;

    Form * getCopy() const ;
} ;


}

#endif // GEL_BEHAVIOUR
