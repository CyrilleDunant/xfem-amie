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

PlasticStrain::PlasticStrain() : previousImposedStrain(0.,3), imposedStrain(0.,3)
{
	getState(true).resize(1, 0.);
	getPreviousState().resize(1, 0.);
	isNull = false ;
	es = NULL ;
	v.push_back(XI);
	v.push_back(ETA);
	param = NULL ;
	plasticVariable = 0 ;
	c_psi = 0.05 ;
	eps_f = 0.0057 ;
}

double PlasticStrain::plasticFlowPotential(const Matrix &m) const
{
	Matrix s(m-identity(m.numCols())*trace(m)/(double)m.numCols()) ;
	return c_psi * trace(m) + sqrt(0.5*s.squareFroebeniusNorm()) ;
}

// void PlasticStrain::step(ElementState & s)
// {
// 	if(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() > POINT_TOLERANCE_2D)
// 	{
// 		if(v.size() == 2 && !previousImposedStrain.size())
// 			previousImposedStrain.resize(3, 0.) ;
// 		
// 		if(v.size() == 3 && !previousImposedStrain.size())
// 			previousImposedStrain.resize(6, 0.) ;
// 		
// 		if(v.size() == 2 && !imposedStrain.size())
// 			imposedStrain.resize(3, 0.) ;
// 		
// 		if(v.size() == 3 && !imposedStrain.size())
// 			imposedStrain.resize(6, 0.) ;
// 		
// 		converged = false ;
// 		Vector istr = imposedStrain ;
// 		computeDamageIncrement(s) ;
// 		converged = std::abs(istr-imposedStrain).max() < 1e-4 ;
// 		change = converged ;
// 	}
// }

std::pair<Vector, Vector> PlasticStrain::computeDamageIncrement(ElementState & s)
{
	if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL && v.size() == 2)
	{
		v.push_back(ZETA);
		previousImposedStrain.resize(6, 0.) ;
		imposedStrain.resize(6, 0.) ;
	}
	
	if(!param)
		param = new Matrix(s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;
	if(!es)
		es =&s ;


		if(s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() && s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() )
		{
			
			Matrix stressMatrix(v.size(), v.size()) ;
			Vector stress = s.getParent()->getBehaviour()->getFractureCriterion()->smoothedStress(s, EFFECTIVE_STRESS) ;
			stressMatrix[0][0] = stress[0] ;
			stressMatrix[1][1] = stress[1] ;
			stressMatrix[0][1] = stress[2] ;
			stressMatrix[1][0] = stress[2] ;
			Matrix incrementalStrainMatrix(stressMatrix.numRows(), stressMatrix.numCols()) ;
			
			for(size_t i = 0 ; i < stressMatrix.numRows() ; i++)
			{
				for(size_t j = 0 ; j < stressMatrix.numCols() ; j++)
				{
					Matrix m_p(stressMatrix) ;
					Matrix m_m(stressMatrix) ;
					m_p[i][j] += std::max(std::abs(stressMatrix[i][j]), 1.)*1e-6 ;
					m_m[i][j] -= std::max(std::abs(stressMatrix[i][j]), 1.)*1e-6 ;
					incrementalStrainMatrix[i][j] = (plasticFlowPotential(m_p)-plasticFlowPotential(m_m))/(2e-6*std::max(std::abs(stressMatrix[i][j]), 1.)) ;
				}
			}
// 			incrementalStrainMatrix.print() ;
// 			exit(0) ;
// 			double str = s.getMaximumVonMisesStress() ;
// 			double down_str = str ;
			imposedStrain[0] = incrementalStrainMatrix[0][0] ;
			imposedStrain[1] = incrementalStrainMatrix[1][1] ;
			imposedStrain[2] = incrementalStrainMatrix[1][0] ;
			
			imposedStrain *=  s.getParent()->getBehaviour()->getFractureCriterion()->getFactors().back()*0.01;
		}
		
	return std::make_pair( Vector(0., 1), Vector(1., 1)) ;

}

Matrix PlasticStrain::apply(const Matrix & m) const
{
	if(fractured())
		return m*0. ;
	return m*(1.-getDamage()) ;
}


Matrix PlasticStrain::applyPrevious(const Matrix & m) const
{
	return m ;
}

std::vector<BoundaryCondition * > PlasticStrain::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<BoundaryCondition * > ret ;
	if(!param || fractured())
		return ret ;
// 	Vector f = VirtualMachine().ieval(Gradient(p_i) *( *param *(imposedStrain*state[0]+previousImposedStrain)), gp, Jinv,v) ;
// 	if(state[0] > POINT_TOLERANCE_2D)
// 	{
// 		std::cerr << "strain = "<< (imposedStrain*getState()[0]+previousImposedStrain)[0] << "  " << (imposedStrain*getState()[0]+previousImposedStrain)[1] << std::endl ;
// 		std::cout << "delta = "<< f[0] << "  " << f[1] << std::endl ;
// 		
// 		exit(0) ;
// 	}
	Vector imp = getImposedStress(*p_i.getPoint()) ; //*param*(imposedStrain*getState()[0]+previousImposedStrain) ;
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
		
	return  (Vector)(*param*(imposedStrain*getState()[0]+previousImposedStrain)) ;
}

double PlasticStrain::getDamage() const
{
	return std::max(std::min(1.-exp(-plasticVariable/eps_f), 1.), 0.) ;
}

bool PlasticStrain::fractured() const 
{
// 	if(fraction < 0)
		return false ;
// 	return getDamage() >= thresholdDamageDensity ;
}

void PlasticStrain::postProcess()
{
	if(converged && es && getState()[0] > POINT_TOLERANCE_2D)
	{
		previousImposedStrain += imposedStrain*getState()[0] ;
		std::cout << " score =  " << es->getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() << ", state = " << getState()[0] << ", dstrain = " <<imposedStrain[0]*getState()[0] << "   " <<imposedStrain[1]*getState()[0] << ", strain = "<< previousImposedStrain[0] << "  " << previousImposedStrain[1] << std::endl ;
		plasticVariable += getState()[0]*sqrt(1./3.+2.*c_psi*c_psi) ;
		state[0] = 0;
		imposedStrain = 0 ;
	}

}

PlasticStrain::~PlasticStrain()
{
	delete param ;
}


}
