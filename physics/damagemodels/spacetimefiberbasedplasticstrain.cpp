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
#include "spacetimefiberbasedplasticstrain.h"
#include "../fracturecriteria/fracturecriterion.h"
#include "../viscoelasticity.h"
#include "../../utilities/parser/function_parser.h"
#include "../logarithmic_creep_with_external_parameters.h"

namespace Amie {

SpaceTimeFiberBasedPlasticStrain::SpaceTimeFiberBasedPlasticStrain(double f, double t, double y, bool fix)  : fibreFraction(f), timeTolerance(t), fixed(fix)
{
    thresholdDamageDensity = y ;
    getState(true).resize(1, 0.);
    isNull = false ;
}

std::pair< Vector, Vector > SpaceTimeFiberBasedPlasticStrain::computeDamageIncrement( Amie::ElementState &s)
{
    return std::make_pair(state, Vector(1., 1)) ;
}

void SpaceTimeFiberBasedPlasticStrain::computeDelta(ElementState & s)
{
    delta = 1.-getState()[0] ;
}

Matrix SpaceTimeFiberBasedPlasticStrain::applyViscous(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;

    if(fractured())
        return m*residualStiffnessFraction ;

    return m ;
}

Matrix SpaceTimeFiberBasedPlasticStrain::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;

    if(fractured())
        return m*residualStiffnessFraction ;

    return m ;
}


bool SpaceTimeFiberBasedPlasticStrain::fractured(int direction) const
{
    if(fraction < 0)
        return false ;
    return getState().max() >= thresholdDamageDensity ;
}


SpaceTimeFiberBasedPlasticStrain::~SpaceTimeFiberBasedPlasticStrain()
{
}

void SpaceTimeFiberBasedPlasticStrain::step( ElementState &s , double maxscore)
{
    elementState = &s ;
    if( fraction < 0 )
    {
        plasticStrain.resize(3) ;
        C.resize(3,3) ;
        plasticStrain = 0 ;        
        C = 0 ;
        if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            fraction = s.getParent()->area() ;
        else
            fraction = s.getParent()->volume();
    }

    converged = true ;
    change = false ;
    if(!s.getParent()->getBehaviour()->getFractureCriterion() || maxscore < 0)
    {
        return ;
    }

    double score = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;

    if(!fractured() && score > 0 && (maxscore - score) < timeTolerance)
    {
        state[0] += fibreFraction ;
        change = true ;
        converged = true ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;

        Point p = s.getParent()->getCenter() ;
        p.setT( 1. -2.*maxscore ) ; 
        C = s.getParent()->getBehaviour()->getTensor( p ) ;
        if(dynamic_cast<Viscoelasticity *>( s.getParent()->getBehaviour() )) { C = dynamic_cast<Viscoelasticity *>( s.getParent()->getBehaviour() )->getElasticTensor( p ) ; }
        bool oriented = false ;
        if(fixed && state[0] > 0) { oriented = true ; }

        if(dynamic_cast<LogarithmicCreepWithExternalParameters *>( s.getParent()->getBehaviour() ))
        {
            if(dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).has("plastic_orientation"))
            {
                orientation = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get("plastic_orientation", values ) ;
                oriented = true ;
            }
        }

        if(!oriented)
        {
            Vector stress = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedField(MECHANICAL_STRAIN_FIELD,  s , p.getT()) ;
            orientation = atan2 ( stress[1], stress[0] ) ;
//            std::cout << stress[0] << " " << stress[1] << " " << stress[2] << " " << orientation << std::endl ;
        }

        Vector increment(3) ; 
        increment = 0 ;
        increment[0] = fibreFraction ;
        plasticStrain += Tensor::rotate2ndOrderTensor2D( increment, orientation ) ;
//        std::cout << orientation << " " << plasticStrain[0] << " " << plasticStrain[1] << " " << plasticStrain[2] << std::endl ;
    }

    return ;
}


DamageModel * SpaceTimeFiberBasedPlasticStrain::getCopy() const
{
    SpaceTimeFiberBasedPlasticStrain * dam = new SpaceTimeFiberBasedPlasticStrain(fibreFraction, timeTolerance, thresholdDamageDensity, fixed) ;
    dam->copyEssentialParameters( this ) ;
    return dam ;
}

std::vector<BoundaryCondition * > SpaceTimeFiberBasedPlasticStrain::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;
    if(C.numCols() == 0 || fractured())
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

Vector SpaceTimeFiberBasedPlasticStrain::getImposedStress(const Point & p) const 
{
    if(fractured()) { return Vector(0., plasticStrain.size()) ; }
    return C*plasticStrain ;
}

Vector SpaceTimeFiberBasedPlasticStrain::getImposedStrain(const Point & p) const 
{
    if(fractured()) { return Vector(0., plasticStrain.size()) ; }
    return plasticStrain ;
}

SpaceTimeFiberBasedPlasticDamage::SpaceTimeFiberBasedPlasticDamage(std::string d, std::string reqs, double f, double t, double y, bool fix)  : SpaceTimeFiberBasedPlasticStrain(f,t,y,fix)
{
    getState(true).resize(2, 0.);
    needExternalVariables = false ;
    
    std::vector<std::string> req ;
    if(reqs.length() > 0)
    {
        size_t i = 0 ;
        size_t pos = reqs.find(',') ;
        while( pos != std::string::npos )
        {
            req.push_back( reqs.substr( i, pos-i ) ) ;
            i = ++pos ;
            pos = reqs.find(',', i) ;
        }
       req.push_back( reqs.substr(i, reqs.length() ) ) ;
    }

    std::map<std::string, char> coordDamage ;
    coordDamage["plastic_strain"] = 'x' ;
    if(req.size() > 0)
    {
        needExternalVariables = true ;
        coordDamage[ req[0] ] = 'y' ;
        yCoord = req[0] ;
    }
    if(req.size() > 1)
    {
        needExternalVariables = true ;
        coordDamage[ req[1] ] = 'z' ;
        zCoord = req[1] ;
    }
    if(req.size() > 2)
    {
        needExternalVariables = true ;
        coordDamage[ req[2] ] = 't' ;
        tCoord = req[2] ;
    }
    if(req.size() > 3)
    {
        needExternalVariables = true ;
        coordDamage[ req[3] ] = 'u' ;
        uCoord = req[3] ;
    }
    if(req.size() > 4)
    {
        needExternalVariables = true ;
        coordDamage[ req[4] ] = 'v' ;
        vCoord = req[4] ;
    }
    if(req.size() > 5)
    {
        needExternalVariables = true ;
        coordDamage[ req[5] ] = 'w' ;
        wCoord = req[5] ;
    }
    damage = FunctionParser::getFunction( d, coordDamage ) ;
}

SpaceTimeFiberBasedPlasticDamage::SpaceTimeFiberBasedPlasticDamage( Function d, std::string y, std::string z, std::string t, std::string u, std::string v, std::string w, double f, double tol, double yield, bool fix ) : SpaceTimeFiberBasedPlasticStrain( f, tol, yield, fix ), damage(d), yCoord(y), zCoord(z), tCoord(t), uCoord(u), vCoord(v), wCoord(w)
{
    if(y.length() > 0 || z.length() > 0 || t.length() > 0 || u.length() > 0 || v.length() > 0 || w.length() > 0) { needExternalVariables = true ;}
}


Matrix SpaceTimeFiberBasedPlasticDamage::applyViscous(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;
    if(fractured())
        return m*residualStiffnessFraction ;

    return m*(1.-state[1]) ;
}

Matrix SpaceTimeFiberBasedPlasticDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;
    if(fractured())
        return m*residualStiffnessFraction ;

    return m*(1.-state[1]) ;
}

Vector SpaceTimeFiberBasedPlasticDamage::getImposedStress(const Point & p) const 
{
    if(fractured()) { return Vector(0., plasticStrain.size()) ; }
    return apply(C,p)*plasticStrain ;
}

DamageModel * SpaceTimeFiberBasedPlasticDamage::getCopy() const
{
    SpaceTimeFiberBasedPlasticDamage * dam = new SpaceTimeFiberBasedPlasticDamage(damage, yCoord, zCoord, tCoord, uCoord, vCoord, wCoord, fibreFraction, timeTolerance, thresholdDamageDensity, fixed) ;
    dam->copyEssentialParameters( this ) ;
    return dam ;
}

void SpaceTimeFiberBasedPlasticDamage::step(ElementState & s, double maxscore)  
{
    SpaceTimeFiberBasedPlasticStrain::step(s, maxscore) ;
    
    if(change)
    {
        double ycoor = 0 ;
        double zcoor = 0 ;
        double tcoor = 0 ;
        double ucoor = 0 ;
        double vcoor = 0 ;
        double wcoor = 0 ;

        if( needExternalVariables )
        {
            if( yCoord.length() > 0 )
                ycoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( yCoord, values ) ;
            if( zCoord.length() > 0 )
                zcoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( zCoord, values ) ;
            if( tCoord.length() > 0 )
                tcoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( tCoord, values ) ;
            if( uCoord.length() > 0 )
                ucoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( uCoord, values ) ;
            if( vCoord.length() > 0 )
                vcoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( vCoord, values ) ;
            if( wCoord.length() > 0 )
                wcoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( wCoord, values ) ;
        }

        state[1] = std::min(1., std::max( state[1], VirtualMachine().eval( damage, state[0], ycoor, zcoor, tcoor, ucoor, vcoor, wcoor ) ) ) ;
    }
}




}

