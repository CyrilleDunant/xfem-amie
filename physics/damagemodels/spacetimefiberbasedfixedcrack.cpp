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
    if(state.max() < POINT_TOLERANCE || abs(m.array()).max() < POINT_TOLERANCE )
        return m ;

    Matrix C = m ;
    C = Tensor::rotate4thOrderTensor2D( C, orientation ) ;
    C[0][0] *= 1.-state[0] ;
    C[0][1] *= 1.-(state.min()) ;
    C[0][2] *= 1.-state[0] ;
    C[1][0] *= 1.-(state.min()) ;
    C[2][0] *= 1.-state[0] ;
    C[1][1] *= 1.-state[1] ;
    C[2][1] *= 1.-state[1] ;
    C[1][2] *= 1.-state[1] ;

/*    if(m.array().max() > 1)
    {

    std::cout << "--------------" << std::endl ;

    std::cout << orientation << "\t" << state[0] << "\t" << state[1] << std::endl ;

    m.print() ;

    std::cout << "------" << std::endl ;

    C.print() ;

    std::cout << "------" << std::endl ;

    Tensor::rotate4thOrderTensor2D( C, -orientation ).print() ; 

    std::cout << "------" << std::endl ;
    }*/

    return Tensor::rotate4thOrderTensor2D( C, -orientation ) ;
}



void SpaceTimeFiberBasedFixedCrack::step( ElementState &s , double maxscore)
{
    elementState = &s ;
    converged = true ;

    if( fraction < 0 )
    {
        if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            fraction = s.getParent()->area() ;
        else
            fraction = s.getParent()->volume();

        if( fixedOrientation  )
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
    }

    change = false ;
    if(!s.getParent()->getBehaviour()->getFractureCriterion() || maxscore < 0)
    {
        return ;
    }

    double score = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;

    #pragma omp critical(taratartoto)
    if(!fractured() && score > 0 && (maxscore - score) < timeTolerance)
    {
        change = true ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;

        if( !fixedOrientation )
        {
            Vector strain = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedField(orientationField,  s , ( 1. -2.*maxscore )) ;
            orientation = 0.5*atan2( 0.5*strain[2], strain[0] - strain[1] ) ;
            fixedOrientation = true ;
        }
        
        Vector strain = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedField(orientationField,  s , ( 1. -2.*maxscore )) ;
        Vector rotated = Tensor::rotate2ndOrderTensor2D( strain, -orientation ) ;
//        std::cout << strain[0] << "\t" << strain[1] << "\t" << strain[2] << "\t" << orientation << "\t" << rotated[0] << "\t"  << rotated[1] << "\t" << rotated[2] << std::endl ;
        if( std::abs( rotated[0] ) > std::abs( rotated[1] ) ) 
            state[0] += fibreFraction ;
        else
            state[1] += fibreFraction ;

        for(size_t i = 0 ; i < state.size() ; i++)
        {
            if(state[i] > 1)
                state[i] = 1. ;
        }
    }

    return ;
}

DamageModel * SpaceTimeFiberBasedFixedCrack::getCopy() const
{
    SpaceTimeFiberBasedFixedCrack * dam = new SpaceTimeFiberBasedFixedCrack(fibreFraction, timeTolerance, thresholdDamageDensity, orientationField, fixedOrientation && (fraction < 0) ) ;
    dam->copyEssentialParameters( this ) ;
    return dam ;
}

bool SpaceTimeFiberBasedFixedCrack::fractured(int direction) const
{
    if(fraction < 0)
        return false ;
    return getState().max() >= thresholdDamageDensity ;
}

}

