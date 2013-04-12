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
#include "plasticstrain.h"
#include "../../features/boundarycondition.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu {

PlasticStrain::PlasticStrain() : previousCompressiveImposedStrain(0.,3), previousTensileImposedStrain(0.,3), imposedStrain(0.,3)
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
	c_psi = 0.05 ;
	eps_f = 0.0057; //0.0057 ;
	kappa_0 = 0 ; //3.350e-3 ;  //3.350e-3 ; // ; //3.350e-3 ; //5 up ; 4 down
	es = nullptr ;
	broken = false ;
	factor = 1 ;

}

double PlasticStrain::plasticFlowPotential(const Matrix &m) const
{
	Matrix s(m-identity(m.numCols())*trace(m)/(double)m.numCols()) ;
	return c_psi * trace(m) + sqrt(0.5*s.squareFroebeniusNorm()) ;
}

std::pair<Vector, Vector> PlasticStrain::computeDamageIncrement(ElementState & s)
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
		setConvergenceType(DISSIPATIVE_CENTER);
	}
	
	if(!param)
		param = new Matrix(s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;
	
	
	if( s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() && s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint())
	{
		
		Vector originalIstrain = getImposedStrain(s.getParent()->getCenter()) ;
		Matrix stressMatrix(v.size(), v.size()) ;
		Vector stress(3) ;
		Vector strain(3) ;
// 		s.getField(STRAIN_FIELD, REAL_STRESS_FIELD,es->getParent()->getCenter(),strain,stress, false);
		std::pair<Vector, Vector> ss = s.getParent()->getBehaviour()->getFractureCriterion()->smoothedStressAndStrain(s) ;
		stress = ss.first ;
		strain = ss.second ;
		stressMatrix[0][0] = stress[0] ;
		stressMatrix[1][1] = stress[1] ;
		stressMatrix[0][1] = stress[2] ;
		stressMatrix[1][0] = stress[2] ;
		Matrix incrementalStrainMatrix(stressMatrix.numRows(), stressMatrix.numCols()) ;
		double iftynorm = std::abs(stressMatrix.array()).max() ;
		if(iftynorm < POINT_TOLERANCE_2D)
			iftynorm = 1 ;
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

		
		imposedStrain /=sqrt(imposedStrain[0]*imposedStrain[0]+imposedStrain[1]*imposedStrain[1]+imposedStrain[2]*imposedStrain[2]) ;
		double inftysnorm = std::abs(strain-originalIstrain).max() ;
		imposedStrain*=inftysnorm ;
		inCompression = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(1) ;
		inTension = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0) ;
	}
		
	return std::make_pair( Vector(0., 1), Vector(1., 1)) ;

}

int PlasticStrain::getMode() const 
{ 
// 	if(es
// 		&& es->getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() 
// 		&&(
// 			inCompression != es->getParent()->getBehaviour()->getFractureCriterion()->directionMet(1)
// 		|| inTension != es->getParent()->getBehaviour()->getFractureCriterion()->directionMet(0) 
// 		|| inCompression && compressivePlasticVariable <= kappa_0 && getPlasticity() >= kappa_0-POINT_TOLERANCE_2D
// 		|| inTension && tensilePlasticVariable <= kappa_0 && getPlasticity() >= kappa_0-POINT_TOLERANCE_2D
// 		)
// 	) 
// 	{
// 		return 1 ;
// 	}
	return -1 ;
}

double PlasticStrain::getAngleShift() const
{
// 	if(!es)
		return 0 ;
		
		Matrix stressMatrix(v.size(), v.size()) ;
		Vector stress(3) ;
		Vector strain(3) ;
		es->getField(STRAIN_FIELD, REAL_STRESS_FIELD,es->getParent()->getCenter(),strain,stress, false);
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
		if(std::abs(imposedStrain).max() > POINT_TOLERANCE_2D)
		{
			istrain /= sqrt(istrain[0]*istrain[0]+istrain[1]*istrain[1]+istrain[2]*istrain[2]) ;
			istrain *= .1;nimposed*(1.-getDamage()) ;
		}
		double dp = std::inner_product(&istrain[0],&istrain[3], &imposedStrain[0], double(0.))/(nimposed*nimposed) ;
		double angle = acos(dp) ;
		if(std::abs(dp-1.) < POINT_TOLERANCE_2D)
			angle = 0 ;
		if(nimposed < POINT_TOLERANCE_2D)
			angle = 0 ;
		if (angle < 0 )
			angle += M_PI ;
// 		std::cout << std::endl ;
// 		std::cout << istrain[0] << "  " << istrain[1] << "  "<< istrain[2] << "  " << std::endl ;
// 		std::cout << imposedStrain[0] << "  " << imposedStrain[1] << "  "<< imposedStrain[2] << "  " << std::endl ;
// 		std::cout << angle << std::endl ;
// 		std::cout << std::endl ;
		return angle ;
}

void PlasticStrain::computeDelta(const ElementState & s)
{
	delta = 1 ;
}

Matrix PlasticStrain::apply(const Matrix & m, const Point & p , const IntegrableEntity * e, int g) const
{
	if(fractured())
		return m*1e-6 ;
	return m*(1.-getDamage()) ;
}

std::vector<BoundaryCondition * > PlasticStrain::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<BoundaryCondition * > ret ;
	if(!param || fractured())
		return ret ;

	Vector imp = getImposedStress(*p_i.getPoint()) ; 
	if(v.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, id, imp[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, id, imp[1]));
//		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, id, imp[2]));
	}
	if(v.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[2]));
/*		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[3]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[5]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp, Jinv, id, imp[4]));*/
	}
	return ret ;
}

Vector PlasticStrain::getImposedStress(const Point & p) const
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
		
	return  (Vector)(*param*(1.-getDamage())*getImposedStrain(p)) ;
}

Vector PlasticStrain::getImposedStrain(const Point & p) const
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
		return  (imposedStrain*getState()[0]+previousCompressiveImposedStrain) ;
	
// 	return  imposedStrain*getState()[0]+previousTensileImposedStrain ;
}

double PlasticStrain::getDamage() const
{

	//return std::min(topdamage*state[0]+bottomdamage*(1.-state[0]) + factor, 1.);
	
	double currentPlaticVariable = getPlasticity() ;
	if(currentPlaticVariable >= kappa_0*factor)
	{
// 		return std::max((currentPlaticVariable-kappa_0)/eps_f,0.) ;
		return 1.-exp(-(currentPlaticVariable-kappa_0*factor)/(eps_f)) ;
	}
	return 0 ;
}

double PlasticStrain::getPlasticity() const
{
	Vector istrain = imposedStrain*getState()[0];
	double currentPlaticVariable = sqrt(2./3.)*sqrt(istrain[0]*istrain[0]+istrain[1]*istrain[1]+istrain[2]*istrain[2]) ;
// 	if(inCompression )
		currentPlaticVariable += compressivePlasticVariable ;// sqrt(2./3.)*sqrt(previousCompressiveImposedStrain[0]*previousCompressiveImposedStrain[0]+previousCompressiveImposedStrain[1]*previousCompressiveImposedStrain[1]+previousCompressiveImposedStrain[2]*previousCompressiveImposedStrain[2]) ;
// 	else
// 		currentPlaticVariable += tensilePlasticVariable ;// sqrt(2./3.)*sqrt(previousTensileImposedStrain[0]*previousTensileImposedStrain[0]+previousTensileImposedStrain[1]*previousTensileImposedStrain[1]+previousTensileImposedStrain[2]*previousTensileImposedStrain[2]) ;
	return currentPlaticVariable ;
}

bool PlasticStrain::fractured() const 
{
	if(fraction < 0)
		return false ;
	return broken || getDamage() >= thresholdDamageDensity ;
}

void PlasticStrain::postProcess()
{
	if(converged && es && state[0] > POINT_TOLERANCE_2D)
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

PlasticStrain::~PlasticStrain()
{
	delete param ;
}


}
