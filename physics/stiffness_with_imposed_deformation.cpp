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

#include "stiffness_with_imposed_deformation.h"

using namespace Mu ;

StiffnessWithImposedDeformation::StiffnessWithImposedDeformation(const Matrix & rig, Vector imposedDef) : LinearForm(rig, false, false, rig.numRows()/3+1) , imposed(imposedDef)
{
	this->time_d = false ;
} ;

StiffnessWithImposedDeformation::~StiffnessWithImposedDeformation() { } ;

Matrix StiffnessWithImposedDeformation::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	
	
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}
Matrix StiffnessWithImposedDeformation::apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	
	
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v) ;
}

bool StiffnessWithImposedDeformation::fractured() const
{
	return false ;
}

Form * StiffnessWithImposedDeformation::getCopy() const 
{
	return new StiffnessWithImposedDeformation(*this) ;
}

bool StiffnessWithImposedDeformation::hasInducedForces() const 
{
	return true ; 
} 

Vector StiffnessWithImposedDeformation::getForces(const ElementState & s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const 
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	
	
	return VirtualMachine().ieval(Gradient(p_i,true) * (param * imposed), gp, Jinv,v) ;
}

