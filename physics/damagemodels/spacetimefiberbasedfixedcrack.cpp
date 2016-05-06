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

    if(fractured(0) && fractured(1))
        return m*residualStiffnessFraction ;
    else
    {
        Matrix S = inverse3x3Matrix( m ) ;
        S = Tensor::rotate4thOrderTensor2D( S, orientation, -1) ;

        Matrix C = inverse3x3Matrix(S) ;

        Matrix D(3,3) ;
        D[0][0] = std::max(std::sqrt(1.-state[0]), residualStiffnessFraction ) ;
        if(cleavage) { D[2][2] = D[0][0] ; D[1][1] = 1 ;}
        else
        {
            D[1][1] = std::max(std::sqrt(1.-state[1]), residualStiffnessFraction ) ;
            D[2][2] = std::sqrt(D[0][0]*D[1][1]) ;
        }

        C = (D*C)*D ;

        if(fractured())
        {
            C = Tensor::rotate4thOrderTensor2D( C, -orientation, 1) ;
            return C ;
        }

//        return C ;

        S = inverse3x3Matrix(C) ;

        S = Tensor::rotate4thOrderTensor2D( S, -orientation, -1 ) ;

        S = inverse3x3Matrix( S ) ;

        return S ;
    }

    return m*residualStiffnessFraction ;
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

    if(!(fractured(0) && fractured(1)) && score > 0 && (maxscore - score) < timeTolerance)
    {
        change = true ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;

        if(state.max() < POINT_TOLERANCE)
        if( !fixedOrientation )
        {
            Vector strain = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedField(orientationField,  s , ( 1. -2.*maxscore )) ;
            orientation = 0.5*atan2( 2*strain[2], strain[1] - strain[0] ) ;
            fixedOrientation = true ;
            if(cleavage)
            {
                Vector rotated = Tensor::rotate2ndOrderTensor2D( strain, orientation ) ;
                if(rotated[1] > rotated[0]) { orientation += M_PI*0.5 ; }
                rotated = Tensor::rotate2ndOrderTensor2D( strain, orientation ) ;
            }
        }
        
        if(cleavage)
        {
            Vector strain = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedField(orientationField,  s , ( 1. -2.*maxscore )) ;
            Vector rotated = Tensor::rotate2ndOrderTensor2D( strain, orientation ) ;
//            std::cout << maxscore << "\t" << rotated[0] << "\t" << rotated[1] << "\t" << rotated[2] << "\t" << state[0] << std::endl ;
            if(rotated[1] > std::max( rotated[0], rotated[2])) { state[0] = 1. ; }
            else
                state[0] += fibreFraction ;
        }
        else
        {
            Vector strain = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedField(orientationField,  s , ( 1. -2.*maxscore )) ;
            Vector rotated = Tensor::rotate2ndOrderTensor2D( strain, orientation ) ;
            if(rotated[2] > std::max( rotated[0], rotated[1]))
            {
                if(state[0] < 1.)
                    state[0] += fibreFraction*std::sqrt(0.5) ;
                if(state[1] < 1.)
                    state[1] += fibreFraction*std::sqrt(0.5) ;
            }
            else if( rotated[0] >  rotated[1] )
            {
                if(state[0] < 1.)
                    state[0] += fibreFraction ;
                else
                    state[1] += fibreFraction ;
            }
            else
            {
                if(state[1] < 1.)
                    state[1] += fibreFraction ;
                else
                    state[0] += fibreFraction ;
            }
        }

        for(size_t i = 0 ; i < state.size() ; i++)
        {
            if(state[i] > thresholdDamageDensity)
                state[i] = 1. ;
        }

    }

    return ;
}

DamageModel * SpaceTimeFiberBasedFixedCrack::getCopy() const
{
    SpaceTimeFiberBasedFixedCrack * dam = new SpaceTimeFiberBasedFixedCrack(fibreFraction, timeTolerance, thresholdDamageDensity, orientationField, fixedOrientation && (fraction < 0), cleavage ) ;
    dam->copyEssentialParameters( this ) ;
    return dam ;
}

bool SpaceTimeFiberBasedFixedCrack::fractured(int direction) const
{
    if(fraction < 0)
        return false ;
    if(cleavage)
    {
        return state[0] >= thresholdDamageDensity ;
    }
    if(direction < 0 || direction >= (int) state.size())
    {
        return this->fractured(0) || this->fractured(1) ;
        return getState().max() >= thresholdDamageDensity ;
    }
    return state[direction] >= thresholdDamageDensity ;
}









}
