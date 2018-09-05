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
#include "spacetimebifurcation.h"
#include "../fracturecriteria/fracturecriterion.h"
#include "../fracturecriteria/maxstrain.h"

namespace Amie {

Matrix SpaceTimeBifurcation::applyViscous(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;

    return m*stiffnessFactor ;
}

Matrix SpaceTimeBifurcation::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;

    return m*stiffnessFactor ;
}

DamageModel * SpaceTimeBifurcation::getCopy() const
{
    SpaceTimeBifurcation * dam = new SpaceTimeBifurcation(stiffnessFactor, timeTolerance) ;
    dam->copyEssentialParameters( this ) ;
    return dam ;
}

std::vector<BoundaryCondition * > SpaceTimeBifurcation::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;
    if(residualStress.size() == 0 || state.max() < POINT_TOLERANCE || fraction < 0)
        return ret ;

    Vector imp = getImposedStress(*p_i.getPoint()) ;
    if(imp.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, id, imp[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, id, imp[1]));
    }
    if(imp.size() == 6)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[2]));
    }
    return ret ;
}

Vector SpaceTimeBifurcation::getImposedStress(const Point & p) const
{
    if(state.max() < POINT_TOLERANCE || fraction < 0) { return Vector(0., residualStrain.size()) ; }
    return residualStress ;
}

Vector SpaceTimeBifurcation::getImposedStrain(const Point & p) const
{
    if(state.max() < POINT_TOLERANCE || fraction < 0) { return Vector(0., residualStrain.size()) ; }
    return residualStrain ;
}

void SpaceTimeBifurcation::step( ElementState &s , double maxscore)
{
    elementState = &s ;
    converged = true ;
    if( fraction < 0 )
    {
        if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            fraction = s.getParent()->area() ;
            residualStrain.resize(3) ;
            residualStress.resize(3) ;
        }
        else
        {
            fraction = s.getParent()->volume();
            residualStrain.resize(6) ;
            residualStress.resize(6) ;
        }
    }

    change = false ;
    if(fractured(0) || !s.getParent()->getBehaviour()->getFractureCriterion() || maxscore < 0)
    {
        return ;
    }

    double score = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;

    if(score > 0 && (maxscore - score) < timeTolerance)
    {
        state[0] = 1. ;
        change = true ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
        Vector strain = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedField( MECHANICAL_STRAIN_FIELD, s, 1.-2.*maxscore) ;
        residualStrain = strain*(1.-1./stiffnessFactor) ;
        Matrix C = dynamic_cast<Viscoelasticity *>(s.getParent()->getBehaviour())->getElasticTensor( Point(0,0,0,1.-2.*maxscore)) ;
        residualStress = (C*stiffnessFactor)*residualStrain ;
/*        std::cout << strain[0] << "\t" << strain[1] << "\t"<< strain[2] << std::endl ;
        std::cout << residualStrain[0] << "\t" << residualStrain[1] << "\t"<< residualStrain[2] << std::endl ;
        std::cout << residualStress[0] << "\t" << residualStress[1] << "\t"<< residualStress[2] << std::endl ;*/
    }

    return ;
}


DamageModel * SpaceTimeBifurcationAndDamage::getCopy() const
{
    SpaceTimeBifurcationAndDamage * dam = new SpaceTimeBifurcationAndDamage( stiffnessFactor, secondary->getCopy(), timeTolerance ) ;
    dam->copyEssentialParameters(this);
    return dam ;
}

bool SpaceTimeBifurcationAndDamage::fractured(int direction) const
{
    if(state[0] > POINT_TOLERANCE)
    {
        if(direction == 0)
            return true ;
        if(direction < 0)
            return secondary->fractured() ;
        return secondary->fractured(direction-1) ;
    }
    else
        return false ;
}

Vector SpaceTimeBifurcationAndDamage::getImposedStrain(const Point &p) const
{
    Vector str = SpaceTimeBifurcation::getImposedStrain(p) ;
    if(secondary->hasInducedForces())
        str += secondary->getImposedStrain(p) ;
    return str ;
}

Vector SpaceTimeBifurcationAndDamage::getImposedStress(const Point &p) const
{
    Vector str = SpaceTimeBifurcation::getImposedStress(p) ;
    if(state[1] > POINT_TOLERANCE)
    {
        Matrix C(str.size(), str.size()) ;
        for(size_t i = 0 ; i < str.size() ; i++)
            C[i][i] = 1 ;
        C = secondary->apply(C,p) ;
        str = C*str ;
    }
    if(secondary->hasInducedForces())
        str += secondary->getImposedStress(p) ;
    return str ;
}

Matrix SpaceTimeBifurcationAndDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;

    Matrix w = secondary->apply(m,p,e,g) ;
    return w*stiffnessFactor ;
}

void SpaceTimeBifurcationAndDamage::step( ElementState &s , double maxscore)
{
    if(!secondary)
    {
        std::cout << "space-time bifurcation: not secondary" << std::endl ;
        exit(0) ;
    }

    elementState = &s ;
    converged = true ;
    if( fraction < 0 )
    {
        if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            fraction = s.getParent()->area() ;
            residualStrain.resize(3) ;
            residualStress.resize(3) ;
        }
        else
        {
            fraction = s.getParent()->volume();
            residualStrain.resize(6) ;
            residualStress.resize(6) ;
        }

        secondary->step(s,-1) ;
        state.resize( secondary->getState().size() +1 ) ;
        state = 0 ;
    }

    change = false ;
    if(fractured() || !s.getParent()->getBehaviour()->getFractureCriterion() || maxscore < 0)
    {
        return ;
    }

    double score = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;

    if(score > 0 && (maxscore - score) < timeTolerance)
    {
        change = true ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
        if( s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0, maxscore) && state[0] < POINT_TOLERANCE)
        {
            state[0] = 1. ;
            Vector strain = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedField( MECHANICAL_STRAIN_FIELD, s, 1.-2.*maxscore) ;
            residualStrain = strain*(1.-1./stiffnessFactor) ;
            Matrix C = dynamic_cast<Viscoelasticity *>(s.getParent()->getBehaviour())->getElasticTensor( Point(0,0,0,1.-2.*maxscore)) ;
            residualStress = (C*stiffnessFactor)*residualStrain ;
        }
        else
        {
            secondary->step(s, maxscore) ;
            for(size_t i = 1 ; i < state.size() ; i++)
            {
                state[i] = secondary->getState()[i-1] ;
            }

        }
    }

    return ;
}


}
