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
    tensilePlasticVariable = 0 ;
    inCompression = false ;
    inTension = false ;
    newtonIteration = false ;
    es = nullptr ;
    broken = false ;
    factor = 0 ;

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
        setConvergenceType(CONSERVATIVE);
    }

    if(!param)
        param = new Matrix(s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;


    if( s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() && s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint())
    {

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
        Matrix m_m(stressMatrix) ;
        double delta = 1e-6*std::abs(plasticFlowPotential(stressMatrix)) ;
        for(size_t i = 0 ; i < stressMatrix.numRows() ; i++)
        {
            for(size_t j = 0 ; j < stressMatrix.numCols() ; j++)
            {

                m_p[i][j] += delta ;
                m_m[i][j] -= delta ;
                incrementalStrainMatrix[i][j] = (plasticFlowPotential(m_p)-plasticFlowPotential(m_m))/(2.*delta) ;
                m_p[i][j] = stressMatrix[i][j] ;
                m_m[i][j] = stressMatrix[i][j] ;
            }
        }
        imposedStrain[0] = incrementalStrainMatrix[0][0] ;
        imposedStrain[1] = incrementalStrainMatrix[1][1] ;
        imposedStrain[2] = 0.5*(incrementalStrainMatrix[0][1]+incrementalStrainMatrix[1][0]) ;


        double norm = sqrt((imposedStrain*imposedStrain).sum()) ;
        double onorm = .5*sqrt(((strain-originalIstrain)*(strain-originalIstrain)).sum())/(*param[0][0]) ;
        if(norm > POINT_TOLERANCE && onorm > POINT_TOLERANCE)
        {
            imposedStrain /= norm ;
            imposedStrain *= onorm ;
        }


        inCompression = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(1) ;
        inTension = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0) ;

        double maxfact = 1. ; //std::max(damageDensityTolerance*5e2, std::min(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()*.5, 1.)) ;
        return std::make_pair( Vector(0., 1), Vector(maxfact, 1)) ;
    }

    return std::make_pair( Vector(0., 1), Vector(1., 1)) ;

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
            double delta = 1e-4*iftynorm ;
            m_p[i][j] += delta ;
            m_m[i][j] -= delta ;
            m_p2[i][j] += 2.*delta ;
            m_m2[i][j] -= 2.*delta ;
            incrementalStrainMatrix[i][j] = ( plasticFlowPotential(m_m2)/12. - 2./3.*plasticFlowPotential(m_m) + 2./3.*plasticFlowPotential(m_p) - plasticFlowPotential(m_p2)/12. ) / (4.*delta) ;
            m_p = stressMatrix ;
            m_m = stressMatrix ;
            m_p2 = stressMatrix ;
            m_m2 = stressMatrix ;
        }
    }

    istrain[0] = incrementalStrainMatrix[0][0] ;
    istrain[1] = incrementalStrainMatrix[1][1] ;
    istrain[2] = 0.5*(incrementalStrainMatrix[0][1]+incrementalStrainMatrix[1][0]) ;
    double nimposed = sqrt(std::inner_product(&imposedStrain[0],&imposedStrain[3],&imposedStrain[0],0.)) ;
    if(std::abs(imposedStrain).max() > POINT_TOLERANCE)
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
            Matrix m_m(stressMatrix) ;
            double delta = 1e-6*std::abs(plasticFlowPotential(stressMatrix)) ;
            for(size_t i = 0 ; i < stressMatrix.numRows() ; i++)
            {
                for(size_t j = 0 ; j < stressMatrix.numCols() ; j++)
                {

                    m_p[i][j] += delta ;
                    m_m[i][j] -= delta ;
                    incrementalStrainMatrix[i][j] = (plasticFlowPotential(m_p)-plasticFlowPotential(m_m))/(2.*delta) ;
                    m_p[i][j] = stressMatrix[i][j] ;
                    m_m[i][j] = stressMatrix[i][j] ;
                }
            }
            imposedStrain[0] = incrementalStrainMatrix[0][0] ;
            imposedStrain[1] = incrementalStrainMatrix[1][1] ;
            imposedStrain[2] = 0.5*(incrementalStrainMatrix[0][1]+incrementalStrainMatrix[1][0]) ;


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

            m_p = (stressMatrix) ;
            m_m = (stressMatrix) ;
            delta = 1e-6*std::abs(plasticFlowPotential(stressMatrix)) ;
            for(size_t i = 0 ; i < stressMatrix.numRows() ; i++)
            {
                for(size_t j = 0 ; j < stressMatrix.numCols() ; j++)
                {

                    m_p[i][j] += delta ;
                    m_m[i][j] -= delta ;
                    incrementalStrainMatrix[i][j] = (plasticFlowPotential(m_p)-plasticFlowPotential(m_m))/(2.*delta) ;
                    m_p[i][j] = stressMatrix[i][j] ;
                    m_m[i][j] = stressMatrix[i][j] ;
                }
            }
            imposedStrain[0] = incrementalStrainMatrix[0][0] ;
            imposedStrain[1] = incrementalStrainMatrix[1][1] ;
            imposedStrain[2] = 0.5*(incrementalStrainMatrix[0][1]+incrementalStrainMatrix[1][0]) ;


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
//     if(fractured())
//         return m*1e-6 ;
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

//     if(fractured())
//     {
//         if(v.size() == 2)
//             return Vector(0., 3) ;
//         return Vector(0., 6) ;
//     }
// 	if(inCompression )
    return  imposedStrain*getState()[0]+previousCompressiveImposedStrain ;

// 	return  imposedStrain*getState()[0]+previousTensileImposedStrain ;
}

double PrandtlGrauertPlasticStrain::getDamage() const
{
//      return 0 ;
    //return std::min(topdamage*state[0]+bottomdamage*(1.-state[0]) + factor, 1.);

    double currentPlaticVariable = getPlasticity() ;
    if(currentPlaticVariable >= kappa_0*factor)
    {
// 		return std::max((currentPlaticVariable-kappa_0)/eps_f,0.) ;
        return 1.-exp(-(currentPlaticVariable-kappa_0*factor)/(eps_f*16.)) ;
    }
    return 0 ;
}

double PrandtlGrauertPlasticStrain::getPlasticity() const
{
    Vector istrain = imposedStrain*getState()[0];
    double currentPlaticVariable = sqrt(2./3.)*sqrt(istrain[0]*istrain[0]+istrain[1]*istrain[1]+istrain[2]*istrain[2]) ;
// 	if(inCompression )
    currentPlaticVariable += compressivePlasticVariable ;// sqrt(2./3.)*sqrt(previousCompressiveImposedStrain[0]*previousCompressiveImposedStrain[0]+previousCompressiveImposedStrain[1]*previousCompressiveImposedStrain[1]+previousCompressiveImposedStrain[2]*previousCompressiveImposedStrain[2]) ;
// 	else
// 		currentPlaticVariable += tensilePlasticVariable ;// sqrt(2./3.)*sqrt(previousTensileImposedStrain[0]*previousTensileImposedStrain[0]+previousTensileImposedStrain[1]*previousTensileImposedStrain[1]+previousTensileImposedStrain[2]*previousTensileImposedStrain[2]) ;
    return currentPlaticVariable ;
}

bool PrandtlGrauertPlasticStrain::fractured(int direction) const
{
//     if(fraction < 0)
    return false ;
//     return broken || getDamage() >= thresholdDamageDensity ;
}

void PrandtlGrauertPlasticStrain::postProcess()
{
    if((converged && es && state[0] > 0) ||
            newtonIteration /*&&
      es->getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() < .05*es->getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() &&
      es->getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet()*/
      )
    {

// 		if(inCompression )
// 		{
        previousCompressiveImposedStrain += imposedStrain * getState()[0] ;
        imposedStrain = imposedStrain * getState()[0] ;
        compressivePlasticVariable += sqrt(2./3.) * sqrt( imposedStrain[0]*imposedStrain[0] +
                                      imposedStrain[1]*imposedStrain[1] +
                                      imposedStrain[2]*imposedStrain[2] ) ;

// 		}
// 		else
// 		{
//
// 			previousTensileImposedStrain += imposedStrain *  getState()[0] + dimposedStrain*(1.-getState()[0]);
// 			imposedStrain = imposedStrain *(getState()[0]/*+1e-4*rand()/RAND_MAX*getState()[0]*/) + dimposedStrain* getState()[0]*(1.-getState()[0]);
// 			tensilePlasticVariable += sqrt(2./3.) * sqrt( imposedStrain[0]*imposedStrain[0] +
// 		                                       imposedStrain[1]*imposedStrain[1] +
// 		                                       imposedStrain[2]*imposedStrain[2] ) ;
// 		}
// 		Vector str = imposedStrain ;
// 		es->getAverageField(STRAIN_FIELD, str) ;
// 		if(std::abs(str).max() > 0.03 /*|| compressivePlasticVariable > 5.*eps_f*/)
// 			broken = true ;

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
