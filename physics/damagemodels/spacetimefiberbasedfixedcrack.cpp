//
// C++ Implementation: isotropiclineardamage
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "spacetimefiberbasedfixedcrack.h"
#include "../fracturecriteria/fracturecriterion.h"
#include "../../polynomial/vm_function_extra.h"
#include "../logarithmic_creep_with_external_parameters.h"

namespace Amie {


Matrix SpaceTimeFiberBasedFixedCrack::applyViscous(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    return apply(m, p, e, g) ;
}

Matrix SpaceTimeFiberBasedFixedCrack::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;

    Vector C(3) ;
    C[0] = m[0][0] ;
    C[1] = 0 ; // m[1][1] ;
    C[2] = 0 ;//(m[1][0]+m[0][1])*0.5 ;

    Vector Crot = Tensor::rotate2ndOrderTensor2D(C, orientation) ;
    Matrix C1(3,3) ;
    C1[0][0] = Crot[0] ;
    C1[0][1] = Crot[2] ;
    C1[1][0] = Crot[2] ;
    C1[1][1] = Crot[1] ;

    Matrix ret = m-C1*(state.max()) ;

    return ret ;
}



void SpaceTimeFiberBasedFixedCrack::step( ElementState &s , double maxscore)
{
    if( fixedOrientation && fraction < 0 )
    {
        if(dynamic_cast<LogarithmicCreepWithExternalParameters *>( s.getParent()->getBehaviour() ))
        {
            if(dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).has("crack_orientation"))
            {
                orientation = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get("crack_orientation", values ) ;
                fixedOrientation = true ;
            }
        }
    }

    SpaceTimeFiberBasedIsotropicLinearDamage::step(s, maxscore) ;

    if(change)
    {
        if(!fixedOrientation)
        {
            Vector stress = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedField(MECHANICAL_STRAIN_FIELD,  s , ( 1. -2.*maxscore )) ;
            orientation = atan2 ( stress[1], stress[0] ) ;
            fixedOrientation = true ;
        }
    }

    return ;
}

DamageModel * SpaceTimeFiberBasedFixedCrack::getCopy() const
{
    SpaceTimeFiberBasedFixedCrack * dam = new SpaceTimeFiberBasedFixedCrack(fibreFraction, timeTolerance, thresholdDamageDensity, fixedOrientation && (fraction < 0) ) ;
    dam->copyEssentialParameters( this ) ;
    return dam ;
}

}

