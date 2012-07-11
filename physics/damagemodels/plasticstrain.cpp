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

PlasticStrain::PlasticStrain() : previousCompressiveImposedStrain(0.,3), previousTensileImposedStrain(0.,3), imposedStrain(0.,3), dimposedStrain(0.,3)
{
	getState(true).resize(1, 0.);
	getPreviousState().resize(1, 0.);
	isNull = false ;
	v.push_back(XI);
	v.push_back(ETA);
	param = NULL ;
	compressivePlasticVariable = 0 ;
	tensilePlasticVariable = 0 ;
	inCompression = false ;
	inTension = false ;
	c_psi = 0.05 ;
	eps_f = 0.0057; //0.0057 ;
	kappa_0 = 3.350e-3 ;  //3.350e-3 ; // ; //3.350e-3 ; //5 up ; 4 down
	es = NULL ;
	damageDensityTolerance = 1e-5 ;
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
		dimposedStrain.resize(6, 0.) ;
	}
	if(!es)
	{
		es = &s ;
		setConvergenceType(DISSIPATIVE_CENTER);
	}
	
	if(!param)
		param = new Matrix(s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;
	
	
	if( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() 
		&& s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() )
	{

		Matrix stressMatrix(v.size(), v.size()) ;
		Vector stress(3) ;
		Vector strain(3) ;
		s.getField(STRAIN_FIELD, REAL_STRESS_FIELD,es->getParent()->getCenter(),strain,stress, false);
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
				incrementalStrainMatrix[i][j] = ( plasticFlowPotential(m_m2)*.25 - 2.*plasticFlowPotential(m_m) + 2.*plasticFlowPotential(m_p) - plasticFlowPotential(m_p2)*.25 ) / (12.*delta) ;
				m_p = stressMatrix ;
				m_m = stressMatrix ;
				m_p2 = stressMatrix ;
				m_m2 = stressMatrix ;
			}
		}
		imposedStrain[0] = incrementalStrainMatrix[0][0] ;
		imposedStrain[1] = incrementalStrainMatrix[1][1] ;
		imposedStrain[2] = 0.5*(incrementalStrainMatrix[0][1]+incrementalStrainMatrix[1][0]) ;
		if(std::abs(imposedStrain).max() > POINT_TOLERANCE_2D)
		{
			imposedStrain /= sqrt(imposedStrain[0]*imposedStrain[0]+imposedStrain[1]*imposedStrain[1]+imposedStrain[2]*imposedStrain[2]) ;
			double snorm =  sqrt(strain[0]*strain[0]+strain[1]*strain[1]+strain[2]*strain[2]) ;
			if(snorm > POINT_TOLERANCE_2D && s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint())
				imposedStrain *= snorm*(1.-getDamage())*.5 ;
		}
// 		imposedStrain[2] = 0 ;
// 		std::cout << imposedStrain[0] << "  " << imposedStrain[1] << "  "<< imposedStrain[2] << "  " << std::endl ;
		inCompression = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(1) ;
		inTension = s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0) ;
	}
		
	return std::make_pair( Vector(0., 1), Vector(1., 1)) ;

}

int PlasticStrain::getMode() const 
{ 
	if(es
		&& es->getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() 
		&&(
			inCompression != es->getParent()->getBehaviour()->getFractureCriterion()->directionMet(1)
		|| inTension != es->getParent()->getBehaviour()->getFractureCriterion()->directionMet(0) 
		|| inCompression && compressivePlasticVariable <= kappa_0 && getPlasticity() > kappa_0
		|| inTension && tensilePlasticVariable <= kappa_0 && getPlasticity() > kappa_0
		)
	) 
	{
		return 1 ;
	}
	return -1 ;
}

double PlasticStrain::getAngleShift() const
{
	if(!es)
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
			istrain *= nimposed ;
		}
		
		double angle = acos(std::inner_product(&istrain[0],&istrain[3], &imposedStrain[0], double(0.))/(nimposed*nimposed)) ;
		if(nimposed < POINT_TOLERANCE_2D)
			angle = M_PI ;
		if (angle < 0 )
			angle += M_PI ;

		return angle ;
}

void PlasticStrain::computeDelta(const ElementState & s)
{
	delta = 1 ;
}

Matrix PlasticStrain::apply(const Matrix & m) const
{
	if(fractured())
		return m*0. ;
	return m*(1.-getDamage()) ;
}


Matrix PlasticStrain::applyPrevious(const Matrix & m) const
{
	if(fractured())
		return m*0. ;
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
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()), id, imp[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()), id, imp[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementarySurface *>(s.getParent()), id, imp[2]));
	}
	if(v.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()), id, imp[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, imp[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, imp[2]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, imp[3]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, imp[5]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, imp[4]));
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
	if(inCompression)
		return  imposedStrain*getState()[0]+dimposedStrain*(1.-getState()[0])+previousCompressiveImposedStrain ;
	
	return  imposedStrain*getState()[0]+dimposedStrain*(1.-getState()[0])+previousTensileImposedStrain ;
}

double PlasticStrain::getDamage() const
{
	double currentPlaticVariable = getPlasticity() ;

	if(currentPlaticVariable >= kappa_0)
	{
		return 1.-exp(-(currentPlaticVariable-kappa_0)/(eps_f)) ;
	}
	return 0 ;
}

double PlasticStrain::getPlasticity() const
{
	Vector istrain = imposedStrain*getState()[0]+dimposedStrain*(1.-getState()[0]) ;
	double currentPlaticVariable = sqrt(2./3.)*sqrt(istrain[0]*istrain[0]+istrain[1]*istrain[1]+istrain[2]*istrain[2]) ;
	if(inCompression)
		currentPlaticVariable += compressivePlasticVariable ;
	else
		currentPlaticVariable += tensilePlasticVariable ;
	return currentPlaticVariable ;
}

bool PlasticStrain::fractured() const 
{
	if(fraction < 0)
		return false ;
	return getDamage() >= thresholdDamageDensity ;
}

void PlasticStrain::postProcess()
{
	if(converged && es)
	{
		
		if(inCompression)
		{
			previousCompressiveImposedStrain += imposedStrain * getState()[0] + dimposedStrain* (1.-getState()[0]);
			imposedStrain = imposedStrain * getState()[0] + dimposedStrain*(1.-getState()[0]);
			compressivePlasticVariable += sqrt(2./3.) * sqrt( imposedStrain[0]*imposedStrain[0] + 
		                                       imposedStrain[1]*imposedStrain[1] + 
		                                       imposedStrain[2]*imposedStrain[2] )  ;
		}
		else
		{
			
			previousTensileImposedStrain += imposedStrain * getState()[0] + dimposedStrain*(1.-getState()[0]);
			imposedStrain = imposedStrain * getState()[0] + dimposedStrain* getState()[0]*(1.-getState()[0]);
			tensilePlasticVariable += sqrt(2./3.) * sqrt( imposedStrain[0]*imposedStrain[0] + 
		                                       imposedStrain[1]*imposedStrain[1] + 
		                                       imposedStrain[2]*imposedStrain[2] )  ;
		}
// 		if(state[0] > POINT_TOLERANCE_2D)
// 		{
// 			std::cout <<std::endl ;
// 			std::cout << previousCompressiveImposedStrain[0] << "   "<< previousCompressiveImposedStrain[1] << "   "<< previousCompressiveImposedStrain[2] << std::endl ;
// 			std::cout << previousTensileImposedStrain[0] << "   "<< previousTensileImposedStrain[1] << "   "<< previousTensileImposedStrain[2] << std::endl ;
// 			std::cout <<std::endl ;
// 		}
		
		state[0] = 0;
		imposedStrain = 0 ;
		dimposedStrain = 0 ;
	}
}

PlasticStrain::~PlasticStrain()
{
	delete param ;
}


}
