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

#include "spatially_distributed_stiffness.h"
#include "physics.h"
#include "mohrcoulomb.h"
#include "maxstrain.h"
#include "ruptureenergy.h"
#include "stiffness_and_fracture.h"
#include "vonmises.h"
#include "mohrcoulomb.h"

using namespace Mu ;

SpatiallyDistributedStiffness::SpatiallyDistributedStiffness(const Matrix & rig, const Matrix & pore, double l,double ca, double cb) : LinearForm(rig, true, true, rig.numRows()/3+1), variability(.2), pore(pore)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);

	length = l ;
	distance = l ;
	criteriona = ca ;
	criterionb = cb ;
} ;

SpatiallyDistributedStiffness::~SpatiallyDistributedStiffness() { } ;

Matrix SpatiallyDistributedStiffness::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}

void SpatiallyDistributedStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool SpatiallyDistributedStiffness::fractured() const
{
	return false ;
}

Form * SpatiallyDistributedStiffness::getCopy() const 
{
	double randomVar = (double)rand()/(double)RAND_MAX ;
//	randomVar = 1.*pow(-log(randomVar),.5) ;
	Matrix newTensor (param*(1.-variability)+param*randomVar*variability) ;
	randomVar = (double)rand()/(double)RAND_MAX ;
	double crita = criteriona*(1.-variability)+criteriona*randomVar*variability ;
	double critb = criterionb*(1.-variability)+criterionb*randomVar*variability ;
	newTensor = pore + (newTensor - pore) * distance / length ;
//	if(randomVar > 0.5)
	if(criteriona > 0)
	  return new StiffnessAndFracture(newTensor, new MohrCoulomb(crita,critb)) ;
	return new Stiffness(newTensor) ;	
//	return new Stiffness(pore) ;
//	return new Stiffness(/*pore*(1.-variability)+*/pore/**randomVar*variability*/) ;
}


void SpatiallyDistributedStiffness::getForces(const ElementState & s, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
}

void SpatiallyDistributedStiffness::setDistance(double d)
{
//  std::cout << d << std::endl ;
  distance = d;
}