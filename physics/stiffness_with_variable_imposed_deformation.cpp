//
// C++ Implementation: stiffness_with_imposed_deformation
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_with_variable_imposed_deformation.h"

using namespace Mu ;

StiffnessWithVariableImposedDeformation::StiffnessWithVariableImposedDeformation(const Matrix & rig, Vector imposedDef) : LinearForm(rig, true, false, rig.numRows()/3+1) , imposed(imposedDef)
{	
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);

	this->time_d = false ;
} ;

StiffnessWithVariableImposedDeformation::~StiffnessWithVariableImposedDeformation() { } ;

Matrix StiffnessWithVariableImposedDeformation::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}
void StiffnessWithVariableImposedDeformation::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool StiffnessWithVariableImposedDeformation::fractured() const
{
	return false ;
}

Form * StiffnessWithVariableImposedDeformation::getCopy() const 
{
	return new StiffnessWithVariableImposedDeformation(*this) ;
}

bool StiffnessWithVariableImposedDeformation::hasInducedForces() const 
{
	return true ; 
} 

void StiffnessWithVariableImposedDeformation::step(double timestep, ElementState & currentState)
{
	double uniformRand = (double)rand()/(double)RAND_MAX ;
	imposed[0] += timestep*uniformRand ;
	imposed[1] += timestep*uniformRand ;
}

Vector StiffnessWithVariableImposedDeformation::getImposedStress(const Point & p) const
{
	return (param * imposed) ;
}

void StiffnessWithVariableImposedDeformation::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
	f = VirtualMachine().ieval(Gradient(p_i) * (param * imposed), gp, Jinv,v) ;
}

