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
#include "prandtlgrauertplasticstrain.h"
#include "../../features/boundarycondition.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Amie {

PrandtlGrauertPlasticStrain::PrandtlGrauertPlasticStrain(double c_psi, double eps_f, double kappa_0) : imposedStrain(0.,3), previousCompressiveImposedStrain(0.,3),  previousTensileImposedStrain(0.,3),c_psi(c_psi), eps_f(eps_f), kappa_0(kappa_0)
{
    getState(true).resize(1, 0.);
    isNull = false ;
    v.push_back(XI);
    v.push_back(ETA);
    param = nullptr ;
    compressivePlasticVariable = 0 ;
    damageDensityTolerance = 2e-2 ;
    tensilePlasticVariable = 0 ;
    thresholdDamageDensity = 1.-damageDensityTolerance ;
    forceDeviatoric=false ;
    inCompression = false ;
    inTension = false ;
    newtonIteration = false ;
    es = nullptr ;
    broken = false ;
    factor = 1. ;

}

double PrandtlGrauertPlasticStrain::plasticFlowPotential(const Matrix &m) const
{
    Matrix s(m-identity(m.numCols())*trace(m)/(double)m.numCols()) ;
    return c_psi * trace(m) + sqrt(0.5*s.squareFroebeniusNorm()) ;
}

std::pair<Vector, Vector> PrandtlGrauertPlasticStrain::computeDamageIncrement(ElementState & s)
{
    if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL && v.size() == 2)
    {
        v.push_back(ZETA);
        previousCompressiveImposedStrain.resize(6, 0.) ;
        previousTensileImposedStrain.resize(6, 0.) ;
        imposedStrain.resize(6, 0.) ;
    }
    if(!es)
    {
        es = &s ;
        setConvergenceType(DISSIPATIVE);
    }

    if(!param)
        param = new Matrix(s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;

//     if( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint())
//         previousCompressiveImposedStrain *= .99 ;

    if( s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() && s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint())
    {
        imposedStrain = 0 ;

        Vector originalIstrain = getImposedStrain(s.getParent()->getCenter()) ;
        Matrix stressMatrix(v.size(), v.size()) ;
        Vector stress(3) ;
        Vector strain(3) ;
// 	s.getField(MECHANICAL_STRAIN_FIELD, EFFECTIVE_STRESS_FIELD,es->getParent()->getCenter(),strain,stress, false);
        state[0] = 0 ;
        s.strainAtGaussPointsSet = false ;
        s.stressAtGaussPointsSet = false ;
        std::pair<Vector, Vector> ss = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedFields( MECHANICAL_STRAIN_FIELD, REAL_STRESS_FIELD, s) ;
        stress = ss.second ;
        strain = ss.first ;
        stressMatrix[0][0] = stress[0] ;
        stressMatrix[1][1] = stress[1] ;
        stressMatrix[0][1] = .5*stress[2] ;
        stressMatrix[1][0] = .5*stress[2] ;
        Matrix incrementalStrainMatrix(stressMatrix.numRows(), stressMatrix.numCols()) ;
        double delta = std::max(1e-12*std::abs(plasticFlowPotential(stressMatrix)), 1e-12) ;
        Matrix m_p(stressMatrix) ;
        Matrix m_m(stressMatrix) ;
        Matrix m_p2(stressMatrix) ;
        Matrix m_m2(stressMatrix) ;
        for(size_t i = 0 ; i < stressMatrix.numRows() ; i++)
        {
            for(size_t j = 0 ; j < stressMatrix.numCols() ; j++)
            {
                m_p[i][j] += delta ;
                m_m[i][j] -= delta ;
                m_p2[i][j] += 2.*delta ;
                m_m2[i][j] -= 2.*delta ;
                incrementalStrainMatrix[i][j] = ( plasticFlowPotential(m_m2)/12. - 2./3.*plasticFlowPotential(m_m) + 2./3.*plasticFlowPotential(m_p) - plasticFlowPotential(m_p2)/12. ) / (4.*delta) ;
                m_p[i][j] = stressMatrix[i][j] ;
                m_m[i][j] = stressMatrix[i][j] ;
                m_p2[i][j] = stressMatrix[i][j] ;
                m_m2[i][j] = stressMatrix[i][j] ;
            }
        }
        
        imposedStrain[0] = incrementalStrainMatrix[0][0] ;
        imposedStrain[1] = incrementalStrainMatrix[1][1] ;
        imposedStrain[2] = (incrementalStrainMatrix[0][1]+incrementalStrainMatrix[1][0]) ;
        state[0] = 0 ;
        s.strainAtGaussPointsSet = false ;
        s.stressAtGaussPointsSet = false ;

        double norm = sqrt((imposedStrain*imposedStrain).sum()) ;
        double onorm =  0.25*sqrt((strain*strain).sum()) ;
        if(norm > POINT_TOLERANCE*POINT_TOLERANCE && onorm > POINT_TOLERANCE*POINT_TOLERANCE)
        {
            imposedStrain /= norm ;
            imposedStrain *= onorm ;
        }
        else
            imposedStrain = previousCompressiveImposedStrain*1e-2 ;
        
        state[0] = 1e-3 ;
        double dscore = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
        if(dscore > s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState())
            imposedStrain *=-1 ;
        
        s.strainAtGaussPointsSet = false ;
        s.stressAtGaussPointsSet = false ;
        ss = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedFields( MECHANICAL_STRAIN_FIELD, REAL_STRESS_FIELD, s) ;
        stress = ss.second ;
        strain = ss.first ;
        stressMatrix[0][0] = stress[0] ;
        stressMatrix[1][1] = stress[1] ;
        stressMatrix[0][1] = .5*stress[2] ;
        stressMatrix[1][0] = .5*stress[2] ;
        m_p = stressMatrix ;
        m_m = stressMatrix ;
        m_p2 = stressMatrix ;
        m_m2 = stressMatrix ;
        for(size_t i = 0 ; i < stressMatrix.numRows() ; i++)
        {
            for(size_t j = 0 ; j < stressMatrix.numCols() ; j++)
            {
                m_p[i][j] += delta ;
                m_m[i][j] -= delta ;
                m_p2[i][j] += 2.*delta ;
                m_m2[i][j] -= 2.*delta ;
                incrementalStrainMatrix[i][j] = ( plasticFlowPotential(m_m2)/12. - 2./3.*plasticFlowPotential(m_m) + 2./3.*plasticFlowPotential(m_p) - plasticFlowPotential(m_p2)/12. ) / (4.*delta) ;
                m_p[i][j] = stressMatrix[i][j] ;
                m_m[i][j] = stressMatrix[i][j] ;
                m_p2[i][j] = stressMatrix[i][j] ;
                m_m2[i][j] = stressMatrix[i][j] ;
            }
        }
        imposedStrain[0] = incrementalStrainMatrix[0][0] ;
        imposedStrain[1] = incrementalStrainMatrix[1][1] ;
        imposedStrain[2] = (incrementalStrainMatrix[0][1]+incrementalStrainMatrix[1][0]) ;
        

        norm = sqrt((imposedStrain*imposedStrain).sum()) ;
        onorm = 0.25*sqrt((strain*strain).sum()) ;
        if(norm > POINT_TOLERANCE*POINT_TOLERANCE && onorm > POINT_TOLERANCE*POINT_TOLERANCE)
        {
            imposedStrain /= norm ;
            imposedStrain *= onorm ;
        } 
        else
            imposedStrain = previousCompressiveImposedStrain*1e-2 ;
        
        s.strainAtGaussPointsSet = false ;
        s.stressAtGaussPointsSet = false ;
        
        dscore = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
        if(dscore > s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState())
            imposedStrain *=-1 ;
//         imposedStrain = previousCompressiveImposedStrain*.005 + imposedStrain*.995 ;
//         dscore = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
//         if( std::abs(dscore) >  std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()))
//             imposedStrain *=-1 ;
        
        state[0] = 0 ;

        inCompression = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(1) ;
        inTension = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0) ;

        double maxfact = std::max(std::min(damageDensityTolerance*1e2, 1.), 
                                  std::min(std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()), 1.)) ;
                                  
                                  
        if(forceDeviatoric)
        {
           double tr = imposedStrain[0] + imposedStrain[1] ;
           imposedStrain[0] -= 0.5*tr ;
           imposedStrain[1] -= 0.5*tr ;
        }
        return std::make_pair( Vector(0., 1), Vector(maxfact, 1)) ;
    }

    return std::make_pair( Vector(.0, 1), Vector(1., 1)) ;

}

int PrandtlGrauertPlasticStrain::getMode() const
{
    return -1 ;
}

double PrandtlGrauertPlasticStrain::getAngleShift() const
{
// 	if(!es)
    return 0 ;

    Matrix stressMatrix(v.size(), v.size()) ;
    Vector stress(3) ;
    Vector strain(3) ;
    es->getField(MECHANICAL_STRAIN_FIELD, REAL_STRESS_FIELD,es->getParent()->getCenter(),strain,stress, false);
    Vector istrain(0.,3) ;
    stressMatrix[0][0] = stress[0] ;
    stressMatrix[1][1] = stress[1] ;
    stressMatrix[0][1] = stress[2] ;
    stressMatrix[1][0] = stress[2] ;
    Matrix incrementalStrainMatrix(stressMatrix.numRows(), stressMatrix.numCols()) ;
    double iftynorm = std::abs(stressMatrix.array()).max() + .1 ;
    Matrix m_p(stressMatrix) ;
    Matrix m_m(stressMatrix) ;
    Matrix m_p2(stressMatrix) ;
    Matrix m_m2(stressMatrix) ;
    for(size_t i = 0 ; i < stressMatrix.numRows() ; i++)
    {
        for(size_t j = 0 ; j < stressMatrix.numCols() ; j++)
        {
            double delta = 1e-6*iftynorm ;
            m_p[i][j] += delta ;
            m_m[i][j] -= delta ;
            m_p2[i][j] += 2.*delta ;
            m_m2[i][j] -= 2.*delta ;
            incrementalStrainMatrix[i][j] = ( plasticFlowPotential(m_m2)/12. - 2./3.*plasticFlowPotential(m_m) + 2./3.*plasticFlowPotential(m_p) - plasticFlowPotential(m_p2)/12. ) / (4.*delta) ;
            m_p[i][j] = stressMatrix[i][j] ;
            m_m[i][j] = stressMatrix[i][j] ;
            m_p2[i][j] = stressMatrix[i][j] ;
            m_m2[i][j] = stressMatrix[i][j] ;
        }
    }

    istrain[0] = incrementalStrainMatrix[0][0] ;
    istrain[1] = incrementalStrainMatrix[1][1] ;
    istrain[2] = 0.5*(incrementalStrainMatrix[0][1]+incrementalStrainMatrix[1][0]) ;
    double nimposed = sqrt(std::inner_product(&imposedStrain[0],&imposedStrain[3],&imposedStrain[0],0.)) ;
    if(std::abs(istrain).max() > POINT_TOLERANCE)
    {
        istrain /= sqrt(istrain[0]*istrain[0]+istrain[1]*istrain[1]+istrain[2]*istrain[2]) ;
        istrain *= .1;
//         nimposed*(1.-getDamage()) ;
    }
    double dp = std::inner_product(&istrain[0],&istrain[3], &imposedStrain[0], double(0.))/(nimposed*nimposed) ;
    double angle = acos(dp) ;
    if(std::abs(dp-1.) < POINT_TOLERANCE)
        angle = 0 ;
    if(nimposed < POINT_TOLERANCE)
        angle = 0 ;
    if (angle < 0 )
        angle += M_PI ;
    return angle ;
}

void PrandtlGrauertPlasticStrain::computeDelta(ElementState & s)
{
    delta = 1 ;
}


void PrandtlGrauertPlasticStrain::step( ElementState &s , double maxscore)
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
                  std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()) < .05*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()) || !s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet()))
        {
            change = false ;
            return ;
        }

        if(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() > .05*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() && s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet())
        {
            if(std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()) < .05*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance())
            {

                change = false ;
                return ;
            }

            change = true ;


            Vector originalIstrain = getImposedStrain(s.getParent()->getCenter()) ;
            Matrix stressMatrix(v.size(), v.size()) ;
            Vector stress(3) ;
            Vector strain(3) ;
// 	s.getField(MECHANICAL_STRAIN_FIELD, EFFECTIVE_STRESS_FIELD,es->getParent()->getCenter(),strain,stress, false);
            std::pair<Vector, Vector> ss = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedFields( MECHANICAL_STRAIN_FIELD, REAL_STRESS_FIELD, s) ;
            stress = ss.first ;
            strain = ss.second ;
            stressMatrix[0][0] = stress[0] ;
            stressMatrix[1][1] = stress[1] ;
            stressMatrix[0][1] = stress[2] ;
            stressMatrix[1][0] = stress[2] ;
            Matrix incrementalStrainMatrix(stressMatrix.numRows(), stressMatrix.numCols()) ;

            Matrix m_p(stressMatrix) ;
            double delta = 1e-6*std::abs(plasticFlowPotential(stressMatrix)) ;
            for(size_t i = 0 ; i < stressMatrix.numRows() ; i++)
            {
                for(size_t j = 0 ; j < stressMatrix.numCols() ; j++)
                {

                    m_p[i][j] += delta ;
                    incrementalStrainMatrix[i][j] = (plasticFlowPotential(m_p)-plasticFlowPotential(stressMatrix))/(delta) ;
                    m_p[i][j] = stressMatrix[i][j] ;
                }
            }
            imposedStrain[0] = incrementalStrainMatrix[0][0] ;
            imposedStrain[1] = incrementalStrainMatrix[1][1] ;
            imposedStrain[2] = .25*(incrementalStrainMatrix[0][1]+incrementalStrainMatrix[1][0]) ;


            double norm = sqrt((imposedStrain*imposedStrain).sum()) ;
            double onorm = factor*std::min(1.,sqrt(((strain-originalIstrain)*(strain-originalIstrain)).sum())) ;
            if(norm > POINT_TOLERANCE && onorm > POINT_TOLERANCE)
            {
                imposedStrain /= norm ;
                imposedStrain *= onorm ;
            }

            state[0] += POINT_TOLERANCE ;
            s.strainAtGaussPointsSet = false ;
            s.stressAtGaussPointsSet = false ;

            ss = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedFields( MECHANICAL_STRAIN_FIELD, REAL_STRESS_FIELD, s) ;
            stress = ss.first ;
            strain = ss.second ;
            stressMatrix[0][0] = stress[0] ;
            stressMatrix[1][1] = stress[1] ;
            stressMatrix[0][1] = stress[2] ;
            stressMatrix[1][0] = stress[2] ;
            incrementalStrainMatrix = Matrix(stressMatrix.numRows(), stressMatrix.numCols()) ;

            m_p = stressMatrix ;
            for(size_t i = 0 ; i < stressMatrix.numRows() ; i++)
            {
                for(size_t j = 0 ; j < stressMatrix.numCols() ; j++)
                {

                    m_p[i][j] += delta ;
                    incrementalStrainMatrix[i][j] = (plasticFlowPotential(m_p)-plasticFlowPotential(stressMatrix))/(delta) ;
                    m_p[i][j] = stressMatrix[i][j] ;
                }
            }
            imposedStrain[0] = incrementalStrainMatrix[0][0] ;
            imposedStrain[1] = incrementalStrainMatrix[1][1] ;
            imposedStrain[2] = 0.25*(incrementalStrainMatrix[0][1]+incrementalStrainMatrix[1][0]) ;


            norm = sqrt((imposedStrain*imposedStrain).sum()) ;
            onorm = factor*std::min(factor,sqrt(((strain-originalIstrain)*(strain-originalIstrain)).sum())) ;
            if(norm > POINT_TOLERANCE && onorm > POINT_TOLERANCE)
            {
                imposedStrain /= norm ;
                imposedStrain *= onorm ;
            }

            s.strainAtGaussPointsSet = false ;
            s.stressAtGaussPointsSet = false ;

            double pstate = state[0]-POINT_TOLERANCE ;
            double schange = 0.001*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ; //0.001*std::min(mindelta, s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()) ;
            state[0] += schange-POINT_TOLERANCE ; //s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()/*+mindelta)*0.5*/ ;

            inCompression = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(1) ;
            inTension = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0) ;
        } else {
            change = false ;
            return ;
        }
    }
}

Matrix PrandtlGrauertPlasticStrain::apply(const Matrix & m, const Point & p , const IntegrableEntity * e, int g) const
{
    if(fractured())
        return m*1e-6 ;
    return m*(1.-getDamage()) ;
}

std::vector<BoundaryCondition * > PrandtlGrauertPlasticStrain::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;
    if(!param || fractured())
        return ret ;

    Vector imp = getImposedStrain(*p_i.getPoint())*apply(*param) ;
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

Vector PrandtlGrauertPlasticStrain::getImposedStress(const Point & p) const
{
    if(v.size() == 2 /*&& !param*/)
        return Vector(0., 3) ;
//     if(v.size() == 3 && !param)
    return Vector(0., 6) ;
//     if(fractured())
//     {
//         if(v.size() == 2)
//             return Vector(0., 3) ;
//         return Vector(0., 6) ;
//     }
//
//     return  (Vector)(*param*(1.-getDamage())*getImposedStrain(p)) ;
}

Vector PrandtlGrauertPlasticStrain::getImposedStrain(const Point & p) const
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
// 	if(inCompression )
    return  imposedStrain*getState()[0]+previousCompressiveImposedStrain ;

// 	return  imposedStrain*getState()[0]+previousTensileImposedStrain ;
}

double PrandtlGrauertPlasticStrain::getDamage() const
{
    if(forceDeviatoric)
        return 0 ;
    //return std::min(topdamage*state[0]+bottomdamage*(1.-state[0]) + factor, 1.);

    double currentPlaticVariable = getPlasticity() ;
//     if(currentPlaticVariable >  0.5)
//         std::cout << currentPlaticVariable << "  " << kappa_0*factor << std::endl ;
    if(currentPlaticVariable >= kappa_0*factor)
    {
// 		return std::max((currentPlaticVariable-kappa_0)/eps_f,0.) ;
        return std::max(1.-exp(-(currentPlaticVariable-kappa_0*factor)/(eps_f)), 0.) ;
    }
    return 0 ;
}

double PrandtlGrauertPlasticStrain::getPlasticity() const
{
    Vector istrain = imposedStrain*getState()[0];
    double currentPlaticVariable = sqrt(2./3.)*sqrt(istrain[0]*istrain[0]+istrain[1]*istrain[1]+istrain[2]*istrain[2]) ;
    currentPlaticVariable += compressivePlasticVariable ;// sqrt(2./3.)*sqrt(previousCompressiveImposedStrain[0]*previousCompressiveImposedStrain[0]
    return currentPlaticVariable ;
}

bool PrandtlGrauertPlasticStrain::fractured(int direction) const
{
    if(fraction < 0)
        return false ;
    return broken || getDamage() >= thresholdDamageDensity ;
}

void PrandtlGrauertPlasticStrain::postProcess()
{
    if((converged && es && state[0] > 0) ||newtonIteration )
    {
        previousCompressiveImposedStrain += imposedStrain * getState()[0];
        imposedStrain = imposedStrain * getState()[0] ;
        compressivePlasticVariable += sqrt(2./3.) * sqrt( imposedStrain[0]*imposedStrain[0] +
                                      imposedStrain[1]*imposedStrain[1] +
                                      imposedStrain[2]*imposedStrain[2] ) ;
        imposedStrain = 0 ;                               

        if(getDamage() >= thresholdDamageDensity)
            broken = true ;
    }
    else if(converged && es->getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet())
    {
        state[0] = 0;
        imposedStrain = 0 ;
    }
}

PrandtlGrauertPlasticStrain::~PrandtlGrauertPlasticStrain()
{
    delete param ;
}

DamageModel * PrandtlGrauertPlasticStrain::getCopy() const
{
    PrandtlGrauertPlasticStrain * ret = new PrandtlGrauertPlasticStrain() ;
    ret->factor = factor ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}


}
