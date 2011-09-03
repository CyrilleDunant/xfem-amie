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
	
	c_psi = 0.05 ;
}

double PlasticStrain::plasticFlowPotential(const Matrix &m) const
{
	Matrix s = m-identity(m.numCols())*trace(m) ;
	return c_psi * trace(m) + sqrt(secondInvariant(s*s)) ;
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
	Vector ret = s.getStrain(s.getParent()->getCenter());
	if(ret.size() > 3 && v.size() == 2)
	{
		v.push_back(ZETA);
		previousImposedStrain.resize(6, 0.) ;
		imposedStrain.resize(6, 0.) ;
	}
	if(!param)
		param = new Matrix(s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;
	if(!es)
		es =&s ;

		double multiplier_top  = 1 ;
		double multiplier_down = 0 ;

		if(s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() && s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet())
		{
			Matrix stressMatrix = s.getStressMatrix(s.getParent()->getCenter()) ;
			Matrix incrementalStrainMatrix(v.size(), v.size()) ;
			
			for(size_t i = 0 ; i < v.size() ; i++)
			{
				for(size_t j = 0 ; j < v.size() ; j++)
				{
					Matrix m_p(stressMatrix) ;
					Matrix m_m(stressMatrix) ;
					m_p[i][j] += 1e-6 ;
					m_m[i][j] -= 1e-6 ;
					incrementalStrainMatrix[i][j] = (plasticFlowPotential(m_p)-plasticFlowPotential(m_m))/2e-6 ;
				}
			}
	// 		strainMatrix.print();
// 			double str = s.getMaximumVonMisesStress() ;
// 			double down_str = str ;
			imposedStrain[0] = incrementalStrainMatrix[0][0] ;
			imposedStrain[1] = incrementalStrainMatrix[1][1] ;
			imposedStrain[2] = incrementalStrainMatrix[1][0] ;
			double inorm = sqrt(imposedStrain[0]*imposedStrain[0]+imposedStrain[1]*imposedStrain[1]+imposedStrain[2]*imposedStrain[2]) ;
			if(inorm > POINT_TOLERANCE_2D)
				imposedStrain /= inorm ;
		}
		
	return std::make_pair( Vector(0., 1), Vector(1., 1)) ;

}

Matrix PlasticStrain::apply(const Matrix & m) const
{
	return m ;
}


Matrix PlasticStrain::applyPrevious(const Matrix & m) const
{
	return m ;
}

std::vector<BoundaryCondition * > PlasticStrain::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<BoundaryCondition * > ret ;
	if(!param || (imposedStrain*state[0]+previousImposedStrain).max() < POINT_TOLERANCE_2D)
		return ret ;
// 	Vector f = VirtualMachine().ieval(Gradient(p_i) *( *param *(imposedStrain*state[0]+previousImposedStrain)), gp, Jinv,v) ;
// 	if(state[0] > POINT_TOLERANCE_2D)
// 	{
// 		std::cerr << "strain = "<< (imposedStrain*getState()[0]+previousImposedStrain)[0] << "  " << (imposedStrain*getState()[0]+previousImposedStrain)[1] << std::endl ;
// 		std::cout << "delta = "<< f[0] << "  " << f[1] << std::endl ;
// 		
// 		exit(0) ;
// 	}
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
		
	return  *param *(imposedStrain*getState()[0]+previousImposedStrain) ;
}

bool PlasticStrain::fractured() const 
{
	if(fraction < 0)
		return false ;
	return getState().max() >= thresholdDamageDensity ;
}

void PlasticStrain::postProcess()
{
	
	if(converged && es && es->getDeltaTime() < POINT_TOLERANCE_2D && getState()[0] > POINT_TOLERANCE_2D)
	{
		previousImposedStrain += imposedStrain*getState()[0] ;
// 		std::cout << "dstrain = " <<imposedStrain[0]*getState()[0] << "   " <<imposedStrain[1]*getState()[0] << ", strain = "<< previousImposedStrain[0] << "  " << previousImposedStrain[1] << std::endl ;
		state[0] = 0;
		imposedStrain = 0 ;
	}

}

PlasticStrain::~PlasticStrain()
{
	delete param ;
}


}
