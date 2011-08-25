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
	Matrix s = m-identity(3)*trace(m) ;
	return c_psi * trace(m) + sqrt(secondInvariant(s*s)) ;
}

std::pair<Vector, Vector> PlasticStrain::computeDamageIncrement(ElementState & s)
{
	Vector ret = s.getStrain(s.getParent()->getCenter());
	if(ret.size() > 3 && v.size() == 2)
		v.push_back(ZETA);
	if(!param)
		param = new Matrix(s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;
	if(!es)
		es =&s ;
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() && s.getParent()->getBehaviour()->getFractureCriterion()->met())
	{
		if(v.size() == 2 && !previousImposedStrain.size())
			previousImposedStrain.resize(3, 0.) ;
		
		if(v.size() == 3 && !previousImposedStrain.size())
			previousImposedStrain.resize(6, 0.) ;
		
		if(v.size() == 2 && !imposedStrain.size())
			imposedStrain.resize(3, 0.) ;
		
		if(v.size() == 3 && !imposedStrain.size())
			imposedStrain.resize(6, 0.) ;
		Matrix stressMatrix = s.getStressMatrix(s.getParent()->getCenter()) ;
		Matrix strainMatrix(v.size(), v.size()) ;
		
		for(size_t i = 0 ; i < v.size() ; i++)
		{
			for(size_t j = 0 ; j < v.size() ; j++)
			{
				Matrix m_p(stressMatrix) ;
				Matrix m_m(stressMatrix) ;
				m_p[i][j] += 1e-5 ;
				m_m[i][j] -= 1e-5 ;
				strainMatrix[i][j] = (plasticFlowPotential(m_p)-plasticFlowPotential(m_m))/2e-5 ;
			}
		}
		strainMatrix.print();
		imposedStrain[0] = strainMatrix[0][0] ;
		imposedStrain[1] = strainMatrix[1][1] ;
		imposedStrain[2] = 0.5*strainMatrix[1][0] ;
		
		
		imposedStrain *= 0.01 ;
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
	double m = 0 ;
	
	std::vector<BoundaryCondition * > ret ;
	if(!param || !imposedStrain.size())
		return ret ;
	Vector f = VirtualMachine().ieval(Gradient(p_i) *( *param *(imposedStrain+previousImposedStrain)), gp, Jinv,v) ;
// 	if(state[0] > POINT_TOLERANCE_2D)
// 	{
// 		std::cout << "strain = "<< (imposedStrain*state[0]+previousImposedStrain)[0] << "  " << (imposedStrain*state[0]+previousImposedStrain)[1] << std::endl ;
// 		std::cout << "delta = "<< f[0] << "  " << f[1] << std::endl ;
// 	}

	if(f.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[1]));
	}
	if(f.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[2]));
	}
	return ret ;
}

Vector PlasticStrain::getImposedStress(const Point & p) const
{
	if(v.size() == 2 && !param)
		return Vector(0., 3) ;
	if(v.size() == 3 && !param)
		return Vector(0., 6) ;
		
	return  *param *(imposedStrain+previousImposedStrain) ;
}

bool PlasticStrain::fractured() const 
{
	if(fraction < 0)
		return false ;
	return state.max() >= 1. ;
}

void PlasticStrain::postProcess()
{
	
	if(converged)
	{
		if(v.size() == 2 && lastStress.size() == 0)
			lastStress.resize(3, 0.);
		if(v.size() == 3 && lastStress.size() == 0)
			lastStress.resize(6, 0.);
		
		lastStress = es->getStress(es->getParent()->getCenter()) ;
	}
	
	if(converged && state[0] > POINT_TOLERANCE_2D)
	{

		previousImposedStrain += imposedStrain ;
		
		state[0] = 0;
		imposedStrain = 0 ;
		 

		wasBroken = false ;
	}

}

PlasticStrain::~PlasticStrain()
{
	delete param ;
}


}
