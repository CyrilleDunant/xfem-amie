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

PrandtlReussPlasticStrain::PrandtlReussPlasticStrain() : imposedStrain(0.,3), previousImposedStrain(0.,3)
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
    factor = 1 ;

}


std::pair<Vector, Vector> PrandtlReussPlasticStrain::computeDamageIncrement(ElementState & s)
{
    if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL && v.size() == 2)
    {
        v.push_back(ZETA);
        previousImposedStrain.resize(6, 0.) ;
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
	Vector stress(previousImposedStrain.size()) ;
	Vector strain(previousImposedStrain.size()) ;
//        s.getField( REAL_STRESS_FIELD, STRAIN_FIELD,Point(), stress, strain, true) ;
	std::pair<Vector, Vector> ss = s.getParent()->getBehaviour()->getFractureCriterion()->getSmoothedFields(REAL_STRESS_FIELD, MECHANICAL_STRAIN_FIELD, s) ;
	stress = ss.first ;
	strain = ss.second ;
	double tr = stress.size()==3?(stress[0]+stress[1]):(stress[0]+stress[1]+stress[2]) ;
	if(stress.size()==3)
	{
	  for(size_t i = 0 ; i < 2 ; i++)
	    stress[i] -= tr*.5 ;
	}
	else
	{
	  for(size_t i = 0 ; i < 3 ; i++)
	    stress[i] -= tr*.33333333333333333 ;
	}
       
//        Vector pstresses = toPrincipal(stress, SINGLE_OFF_DIAGONAL_VALUES) ;

	imposedStrain = stress ;
//         imposedStrain[0] = pstresses[0] ;
//         imposedStrain[1] = pstresses[1] ;
	if(stress.size()==3)
	  imposedStrain[2] *= 2. ;
	else
	{
	  imposedStrain[3] *= 2. ;
	  imposedStrain[4] *= 2. ;
	  imposedStrain[5] *= 2. ;
	}


        double norm = sqrt((imposedStrain*imposedStrain).sum()) ;
        imposedStrain /= norm ;
	imposedStrain *= sqrt(((strain-originalIstrain)*(strain-originalIstrain)).sum()) ;
	
	double mini = -1 ;
	double mins = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
	for(double i = -1 ; i < 1 ; i+=0.1)
	{
	  state[0] = i ;
// 	  std::cout << i << "  "<< s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) << "  " << s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() << std::endl ;
	  if(s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) < mins)
	  {
	    mini = i ;
	    mins = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
	  }
	}
// 	exit(0) ;
	if(mini < 0)
	  imposedStrain *= -1; 
	state[0] = 0 ;
	
        inCompression = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(1) ;
        inTension = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0) ;
	 return std::make_pair( Vector(0., 1), Vector(1., 1)) ;
    }

    return std::make_pair( Vector(0., 1), Vector(1., 1)) ;

}

int PrandtlReussPlasticStrain::getMode() const
{

    return -1 ;
}

double PrandtlReussPlasticStrain::getAngleShift() const
{
    return 0 ;

}

void PrandtlReussPlasticStrain::computeDelta(ElementState & s)
{
    delta = 1 ;
}

Matrix PrandtlReussPlasticStrain::apply(const Matrix & m, const Point & p , const IntegrableEntity * e, int g) const
{
    if(fractured())
        return m*1e-6 ;
    return m*(1.-getDamage()) ;
}

std::vector<BoundaryCondition * > PrandtlReussPlasticStrain::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;
    if(!param || fractured())
        return ret ;
    
    Vector imp = getImposedStrain(*p_i.getPoint())*(*param) ;
    if(v.size() == 2)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, id, imp[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, id, imp[1]));
    }
    if(v.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[2]));

    }
    return ret ;
}

Vector PrandtlReussPlasticStrain::getImposedStress(const Point & p) const
{
    if(v.size() == 2 )
        return Vector(0., 3) ;
    if(v.size() == 3 )
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
    return  (imposedStrain*getState()[0]+previousImposedStrain) ;

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
    if(converged && es && state[0] > POINT_TOLERANCE)
    {

        previousImposedStrain += imposedStrain * getState()[0] ;
        imposedStrain = imposedStrain * getState()[0] ;
        plasticVariable += sqrt(2./3.) * sqrt( imposedStrain[0]*imposedStrain[0] +
                                      imposedStrain[1]*imposedStrain[1] +
                                      imposedStrain[2]*imposedStrain[2] ) ;

        state[0] = 0;
        imposedStrain = 0 ;
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
