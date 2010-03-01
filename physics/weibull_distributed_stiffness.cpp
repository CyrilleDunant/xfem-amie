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
#include "maxstrain.h"
#include "ruptureenergy.h"

using namespace Mu ;

WeibullDistributedStiffness::WeibullDistributedStiffness(const Matrix & rig, double cri) : LinearForm(rig, true, true, rig.numRows()/3+1), variability(.1)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);

	criterion = cri ;
} ;

WeibullDistributedStiffness::~WeibullDistributedStiffness() { } ;

Matrix WeibullDistributedStiffness::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}

void WeibullDistributedStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool WeibullDistributedStiffness::fractured() const
{
	return false ;
}

Form * WeibullDistributedStiffness::getCopy() const 
{
	double randomVar = (double)rand()/(double)RAND_MAX ;
//	randomVar = 1.*pow(-log(randomVar),.5) ;
	Matrix newTensor (param*(1.-variability)+param*randomVar*variability*2) ;
	return new StiffnessAndFracture(newTensor, 
		new MohrCoulomb(criterion*(1.-variability)+criterion*randomVar*variability,
		 -8.*(criterion*(1.-variability)+criterion*randomVar*variability))) ;
}


void WeibullDistributedStiffness::getForces(const ElementState & s, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
}

