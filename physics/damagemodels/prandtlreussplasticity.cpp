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
#include "prandtlreussplasticity.h"
#include "../../features/boundarycondition.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Amie {

PrandtlReussPlasticStrain::PrandtlReussPlasticStrain() : imposedStrain(0.,3), previousImposedStrain(0.,3), imposedStrainAccumulator(0.,3)
{
    getState(true).resize(1, 0.);
    isNull = false ;
    v.push_back(XI);
    v.push_back(ETA);
    param = nullptr ;
    plasticVariable = 0 ;
    inCompression = false ;
    inTension = false ;
    es = nullptr ;
    broken = false ;
    factor = .5 ;
    allowBackSearch = false ;
    newtonIteration = false ;

}


std::pair<Vector, Vector> PrandtlReussPlasticStrain::computeDamageIncrement(ElementState & s)
{
    if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL && v.size() == 2)
    {
        v.push_back(ZETA);
        previousImposedStrain.resize(6, 0.) ;
        imposedStrainAccumulator.resize(6, 0.) ;
        imposedStrain.resize(6, 0.) ;
    }

    if(!es)
    {
        es = &s ;
//         setConvergenceType(DISSIPATIVE);
    }

    if(!param)
        param = new Matrix(s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;


    if( s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() /*&& s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint()*/)
    {
        double initialState = state[0] ;
        Vector originalIstrain = getImposedStrain(s.getParent()->getCenter()) ;
        Vector stress(previousImposedStrain.size()) ;
        Vector strain(previousImposedStrain.size()) ;
        std::pair<Vector,Vector> str = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedFields( REAL_STRESS_FIELD,MECHANICAL_STRAIN_FIELD, s )  ;
        stress = str.first ;
        strain = str.second ;
//         s.getField( REAL_STRESS_FIELD, MECHANICAL_STRAIN_FIELD,Point(), stress, strain, true) ;

        double tr = (stress.size()==3)?(stress[0]+stress[1]):(stress[0]+stress[1]+stress[2]) ;
        double trstra = (strain.size()==3)?(strain[0]+strain[1]):(strain[0]+strain[1]+strain[2]) ;
        if(stress.size()==3)
        {
            for(size_t i = 0 ; i < 2 ; i++)
            {
                stress[i] -= tr/2. ;
                strain[i] -= trstra/2. ;
            }

            stress[2] *= .5 ;
        }
        else
        {
            for(size_t i = 0 ; i < 3 ; i++)
            {
                stress[i] -= tr/3. ;
                strain[i] -= trstra/3. ;
            }

            stress[3] *= .5 ;
            stress[4] *= .5 ;
            stress[5] *= .5 ;
        }

        imposedStrain = stress ;

        double norm = sqrt((imposedStrain*imposedStrain).sum()) ;
        double onorm = factor*std::abs(strain-originalIstrain).max() ;
        if(norm > POINT_TOLERANCE )
        {
            imposedStrain /= norm ;
        }
        if(onorm > POINT_TOLERANCE)
            imposedStrain *= onorm ;

        //we need to flow in the right direction
        state[0] += POINT_TOLERANCE ;
        s.strainAtGaussPointsSet = false ;
        s.stressAtGaussPointsSet = false ;
//         s.getField( REAL_STRESS_FIELD, MECHANICAL_STRAIN_FIELD,Point(), stress, strain, true) ;
        str = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedFields( REAL_STRESS_FIELD,MECHANICAL_STRAIN_FIELD, s )  ;
        stress = str.first ;
        strain = str.second ;

        tr = (stress.size()==3)?(stress[0]+stress[1]):(stress[0]+stress[1]+stress[2]) ;
        if(stress.size()==3)
        {
            for(size_t i = 0 ; i < 2 ; i++)
                stress[i] -= tr/2. ;

            stress[2] *= .5 ;
        }
        else
        {
            for(size_t i = 0 ; i < 3 ; i++)
                stress[i] -= tr/3. ;

            stress[3] *= .5 ;
            stress[4] *= .5 ;
            stress[5] *= .5 ;
        }

        imposedStrain = stress ;

        norm = sqrt((imposedStrain*imposedStrain).sum()) ;
        onorm = factor*std::abs(strain-originalIstrain).max() ;
        if(norm > POINT_TOLERANCE )
        {
            imposedStrain /= norm ;
        }
        else
        {
            imposedStrain = 0 ;
        }
        if(onorm > POINT_TOLERANCE)
            imposedStrain *= onorm ;
            
// 	std::cout << imposedStrain[0] << "  " << imposedStrain[1] << "  " << imposedStrain[2] << std::endl ;
//
        state[0] = initialState ;

        s.strainAtGaussPointsSet = false ;
        s.stressAtGaussPointsSet = false ;

        inCompression = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(1) ;
        inTension = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0) ;
        double maxfact = std::max(damageDensityTolerance, std::min(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState(), 1.)) ;
        
        for(size_t i = 0 ; i < imposedStrain.size() ; i++)
        {
            if(std::isnan(imposedStrain[i]))
               imposedStrain[i] = 0 ;        
        }
        
        return std::make_pair( Vector(0., 1), Vector(maxfact, 1)) ;
    }

    return std::make_pair( Vector(0., 1), Vector(1., 1)) ;

}

void PrandtlReussPlasticStrain::step( ElementState &s , double maxscore)
{

    if(!newtonIteration)
        DamageModel::step(s, maxscore) ;
    else
    {
        converged = true ;
        std::pair<double, double> delta = s.getParent()->getBehaviour()->getFractureCriterion()->setChange( s , maxscore) ;
        double mindelta = std::abs(delta.first) > std::abs(delta.second) ? delta.first :  delta.second ;

        computeDamageIncrement(s) ;

        if(std::abs(maxscore) < .5*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() &&
                ((s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() &&
                  std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()) < .05*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()) || 
                  !s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet()))
        {
            change = false ;
            return ;
        }

        if(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() > .05*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() && 
           s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet())
        {

            change = true ;

            Vector originalIstrain = getImposedStrain(s.getParent()->getCenter()) ;
            Vector stress(previousImposedStrain.size()) ;
            Vector strain(previousImposedStrain.size()) ;

//         s.getField( REAL_STRESS_FIELD, MECHANICAL_STRAIN_FIELD,Point(), stress, strain, true) ;
            std::pair<Vector, Vector> str = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedFields( REAL_STRESS_FIELD,MECHANICAL_STRAIN_FIELD, s )  ;
            stress = str.first ;
            strain = str.second ;

            double tr = (stress.size()==3)?(stress[0]+stress[1]):(stress[0]+stress[1]+stress[2]) ;
            double trstra = (strain.size()==3)?(strain[0]+strain[1]):(strain[0]+strain[1]+strain[2]) ;
            if(stress.size()==3)
            {
                for(size_t i = 0 ; i < 2 ; i++)
                {
                    stress[i] -= tr/3. ;
                    strain[i] -= trstra/3. ;
                }

                stress[2] *= 0.5 ;
            }
            else
            {
                for(size_t i = 0 ; i < 3 ; i++)
                {
                    stress[i] -= tr/3. ;
                    strain[i] -= trstra/3. ;
                }

                stress[3] *= 0.5 ;
                stress[4] *= 0.5 ;
                stress[5] *= 0.5 ;
            }

            imposedStrain = stress ;

            double norm = sqrt((imposedStrain*imposedStrain).sum()) ;
            double onorm = factor*std::min(1.,sqrt(((strain-originalIstrain)*(strain-originalIstrain)).sum())) ;
            if(norm > POINT_TOLERANCE )
            {
                imposedStrain /= norm ;
                imposedStrain *= onorm ;
            }

            //we need to flow in the right direction
            state[0] += POINT_TOLERANCE ;
            s.strainAtGaussPointsSet = false ;
            s.stressAtGaussPointsSet = false ;
// 	s.getField( REAL_STRESS_FIELD, MECHANICAL_STRAIN_FIELD,Point(), stress, strain, true) ;
            str = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedFields( REAL_STRESS_FIELD,MECHANICAL_STRAIN_FIELD, s )  ;
            stress = str.first ;
            strain = str.second ;

            tr = (stress.size()==3)?(stress[0]+stress[1]):(stress[0]+stress[1]+stress[2]) ;
            if(stress.size()==3)
            {
                for(size_t i = 0 ; i < 2 ; i++)
                    stress[i] -= tr/3. ;

                stress[2] *= 0.5 ;
            }
            else
            {
                for(size_t i = 0 ; i < 3 ; i++)
                    stress[i] -= tr/3. ;

                stress[3] *= 0.5 ;
                stress[4] *= 0.5 ;
                stress[5] *= 0.5 ;
            }

            imposedStrain = stress ;

            norm = sqrt((imposedStrain*imposedStrain).sum()) ;
            onorm = factor*std::min(1.,sqrt(((strain-originalIstrain)*(strain-originalIstrain)).sum())) ;
            if(norm > POINT_TOLERANCE )
            {
                imposedStrain /= norm ;
                imposedStrain *= onorm ;
            }

            state[0] += std::min(mindelta, s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()) - POINT_TOLERANCE; //0.5*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()  - POINT_TOLERANCE ;0.5*; // ;///*+mindelta)*0.5*/ ;  s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()/*+mindelta)*0.5*/ ;

            inCompression = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(1) ;
            inTension = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0) ;
            
            for(size_t i = 0 ; i < imposedStrain.size() ; i++)
            {
                if(std::isnan(imposedStrain[i]))
                imposedStrain[i] = 0 ;        
            }
        } else {
            change = false ;
            return ;
        }
    }
}

int PrandtlReussPlasticStrain::getMode() const
{
    return -1 ;
}

double PrandtlReussPlasticStrain::getAngleShift() const
{
    if(!es)
        return 0 ;

    Vector stress(previousImposedStrain.size()) ;
    Vector strain(previousImposedStrain.size()) ;
    #pragma omp critical
    es->getField( REAL_STRESS_FIELD, MECHANICAL_STRAIN_FIELD,Point(), stress, strain, true) ;

    double tr = (stress.size()==3)?(stress[0]+stress[1]):(stress[0]+stress[1]+stress[2]) ;
    if(stress.size()==3)
    {
        for(size_t i = 0 ; i < 2 ; i++)
            stress[i] -= tr/2. ;

        stress[2] *= .5 ;
    }
    else
    {
        for(size_t i = 0 ; i < 3 ; i++)
            stress[i] -= tr/3. ;

        stress[3] *= .5 ;
        stress[4] *= .5 ;
        stress[5] *= .5 ;
    }

    return acos((stress*imposedStrain).sum()/sqrt((imposedStrain*imposedStrain).sum()*(stress*stress).sum())) ;

}

void PrandtlReussPlasticStrain::computeDelta(ElementState & s)
{
    delta = 1 ;
}

Matrix PrandtlReussPlasticStrain::apply(const Matrix & m, const Point & p , const IntegrableEntity * e, int g) const
{
    if(fractured())
        return m*1e-6 ;
    return m ;
}

std::vector<BoundaryCondition * > PrandtlReussPlasticStrain::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;
    if(!param || fractured() || hasInducedForces() == false)
        return ret ;

    Vector imp = (*param)*(1.-getDamage())*getImposedStrain(*p_i.getPoint()) ;//getImposedStress(*p_i.getPoint()) ;//
    if(v.size() == 2)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, id, imp[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, id, imp[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, id, imp[2]));
    }
    if(v.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[2]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[3]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[4]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[5]));

    }
    return ret ;
}

Vector PrandtlReussPlasticStrain::getImposedStress(const Point & p) const
{
    if(v.size() == 2 /*&& !param*/)
        return Vector(0., 3) ;
//     if(v.size() == 3 /*&& !param*/)
    return Vector(0., 6) ;

//     if(fractured())
//     {
//         if(v.size() == 2)
//             return Vector(0., 3) ;
//         return Vector(0., 6) ;
//     }
//
//     Vector effImposed = imposedStrain*getState()[0]+previousImposedStrain ;
// //     effImposed[2] = sqrt(effImposed[0]*effImposed[0]) ;
//     return  (Vector)(*param*(1.-getDamage())*effImposed) ;
}

Vector PrandtlReussPlasticStrain::getImposedStrain(const Point & p) const
{
    if(v.size() == 2 && !param)
        return Vector(0., 3) ;
    if(v.size() == 3 && !param)
        return Vector(0., 6) ;

    if(fractured())
    {
        if(v.size() == 2)
            return Vector(0., 3) ;
        return Vector(0., 6) ;
    }
    return  imposedStrain*getState()[0]+previousImposedStrain ;

}

double PrandtlReussPlasticStrain::getDamage() const
{
    return 0 ;
//   double eps_f = 10 ;
//     double currentPlaticVariable = getPlasticity() ;
//     if(currentPlaticVariable >= factor)
//     {
//         return 1.-exp(-(currentPlaticVariable-factor)/(eps_f)) ;
//     }
//
}

double PrandtlReussPlasticStrain::getPlasticity() const
{
    Vector istrain = imposedStrain*getState()[0];
    double currentPlaticVariable = sqrt(2./3.)*sqrt(istrain[0]*istrain[0]+istrain[1]*istrain[1]+istrain[2]*istrain[2]) ;
    currentPlaticVariable += plasticVariable ;
    return currentPlaticVariable ;
}

bool PrandtlReussPlasticStrain::fractured(int direction) const
{
    if(fraction < 0)
        return false ;
    return broken || getDamage() >= thresholdDamageDensity ;
}

void PrandtlReussPlasticStrain::postProcess()
{
    if((converged && es && state[0] > 0) ||
            newtonIteration /*&&
      es->getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() < .05*es->getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() &&
      es->getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet()*/
      )
    {

        Vector delta = imposedStrain * getState()[0] - imposedStrainAccumulator ;
        Vector intermediateSum = previousImposedStrain+delta ;
        imposedStrainAccumulator = (intermediateSum-previousImposedStrain)-delta ;
        previousImposedStrain = intermediateSum ;
        imposedStrain = imposedStrain * getState()[0] ;
        plasticVariable += sqrt(2./3.) * sqrt( imposedStrain[0]*imposedStrain[0] +
                                               imposedStrain[1]*imposedStrain[1] +
                                               imposedStrain[2]*imposedStrain[2] ) ;
//         if(plasticVariable > 1e-2)
//             broken = true ;
// 	factor = std::min(std::max(1e-6, state[0]*2.*.25+factor*.75), 0.5) ;
        state[0] = 0;
        imposedStrain = 0 ;

    }
}

void PrandtlReussPlasticStrain::preProcess( double timeStep, ElementState & currentState )
{
    if( state[0] > POINT_TOLERANCE ) {
        currentState.getParent()->behaviourForcesUpdated = true ;
    }
    else  {
        currentState.getParent()->behaviourForcesUpdated = false ;
    }
}

PrandtlReussPlasticStrain::~PrandtlReussPlasticStrain()
{
    delete param ;
}

DamageModel * PrandtlReussPlasticStrain::getCopy() const
{
    PrandtlReussPlasticStrain * ret = new PrandtlReussPlasticStrain() ;
    ret->factor = factor ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}


}
