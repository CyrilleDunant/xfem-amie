//
// C++ Interface: stiffness_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_and_fracture.h"
#include "../delaunay.h"


using namespace Mu ;


StiffnessAndFracture::StiffnessAndFracture(const Matrix & rig, FractureCriterion * crit) : LinearForm(rig, false, true, rig.numRows()/3+1) 
{
	criterion = crit ;
	frac = false ;
} ;

StiffnessAndFracture::~StiffnessAndFracture() 
{ 
	delete criterion ;
} ;

Matrix StiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36 )
		v.push_back(ZETA);
	
	VirtualMachine vm ;
	return vm.ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}

Matrix StiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v) ;
}

void StiffnessAndFracture::step(double timestep, ElementState * currentState) 
{
// 	std::vector<DelaunayTriangle *> neighbourhood = dynamic_cast<DelaunayTriangle *>(currentState->getParent())->neighbourhood ;
// 	
// 	double maxNeighbourhoodStress = 0 ;
// 	double minNeighbourhoodStress = 0 ;
// 	if(!neighbourhood.empty())
// 	{
// 		for(size_t i = 0 ; i< neighbourhood.size() ; i++)
// 		{
// 			Vector pstress = neighbourhood[i]->getState()->getPrincipalStresses(neighbourhood[i]->getCenter()) ;
// 			double maxStress = pstress.max() ;
// 			double minStress = pstress.min() ;
// 			if(maxStress > maxNeighbourhoodStress)
// 				maxNeighbourhoodStress = maxStress ;
// 			if(minStress < minNeighbourhoodStress)
// 				minNeighbourhoodStress = minStress ;
// 		}
// 	}
// 	
// 	Vector pstress = currentState->getPrincipalStresses(currentState->getParent()->getCenter()) ;
// 	double maxStress = pstress.max();
// 	double minStress = pstress.min();
	if(!frac && criterion->met(currentState) )
	{
// 		if(maxNeighbourhoodStress < maxStress)
// 		{
			frac = true ;
			this->param *= 1e-6 ;
		}
// 	}
/*	
	if(!frac && minStress < -4.*sigmaRupt )
	{
		if(minNeighbourhoodStress > minStress)
		{
			frac = true ;
			this->param *= 1e-6 ;
		}
	}*/
}

bool StiffnessAndFracture::fractured() const
{
	return frac;
}

Form * StiffnessAndFracture::getCopy() const 
{
	return new StiffnessAndFracture(*this) ;
}

Vector StiffnessAndFracture::getForces(const ElementState * s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const 
{
	return Vector(0) ;
}

