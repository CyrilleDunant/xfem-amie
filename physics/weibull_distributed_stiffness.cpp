//
// C++ Implementation: weibull_distributed_stiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "weibull_distributed_stiffness.h"
#include "physics.h"
#include "mohrcoulomb.h"
#include "ruptureenergy.h"

using namespace Mu ;

WeibullDistributedStiffness::WeibullDistributedStiffness(const Matrix & rig, double cri) : LinearForm(rig, true, true, rig.numRows()/3+1) 
{
	criterion = cri ;
} ;

WeibullDistributedStiffness::~WeibullDistributedStiffness() { } ;

Matrix WeibullDistributedStiffness::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	
	
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}
Matrix WeibullDistributedStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	
	
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v) ;
}

bool WeibullDistributedStiffness::fractured() const
{
	return false ;
}

Form * WeibullDistributedStiffness::getCopy() const 
{
	double randomVar = (double)random()/(double)RAND_MAX ;
	randomVar = 1.1284*pow(-log(randomVar),1./2.) ;
	Matrix newTensor = param*randomVar ;
	return new StiffnessAndFracture(newTensor, new MohrCoulomb(criterion*randomVar, -8.*criterion*randomVar)) ;
}


Vector WeibullDistributedStiffness::getForces(const ElementState & s, const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{
	return Vector(0) ;
}

