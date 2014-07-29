//
// C++ Implementation: weibull_distributed_stiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "weibull_distributed_stiffness.h"
#include "physics_base.h"
#include "stiffness.h"
#include "stiffness_and_fracture.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "fracturecriteria/confinedmohrcoulomb.h"
#include "fracturecriteria/mcft.h"
#include "fracturecriteria/confinedmohrcoulombwithstrain.h"
#include "fracturecriteria/maxstrain.h"
#include "fracturecriteria/ruptureenergy.h"
#include "homogenization/homogenization_base.h"
#include "../utilities/random.h"

using namespace Amie ;

WeibullDistributedStiffness::WeibullDistributedStiffness(double E, double nu, SpaceDimensionality dim, double down, double up, planeType pt ,MirrorState mirroring, double dx, double dy, double dz) : LinearForm(Material::cauchyGreen(std::make_pair(E,nu), true,dim, pt), true, true, dim), variability(.1), up(up), down(down), E(E), nu(nu), dim(dim), mirroring(mirroring), dx(dx), dy(dy), dz(dz)
{
	materialRadius = .001;
		
	v.push_back(XI);
	v.push_back(ETA);
	if(dim == SPACE_THREE_DIMENSIONAL)
		v.push_back(ZETA);
	damageModel = nullptr ;
} ;

WeibullDistributedStiffness::~WeibullDistributedStiffness() { } ;

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
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1 - variability + variability*weib ;
	StiffnessAndFracture * copy = new StiffnessAndFracture(
								Material::cauchyGreen(std::make_pair(E,nu), true,dim)*factor, 
								new NonLocalMCFT(
										down*factor ,
										E,
										materialRadius, UPPER_BOUND,
										mirroring, dx, dy, dz)
									 ) ;
	if(damageModel)
	{
		delete copy->dfunc ;
		copy->dfunc = damageModel->getCopy() ;
	}
	copy->criterion->setMaterialCharacteristicRadius(materialRadius);
	copy->dfunc->setThresholdDamageDensity(.99);
	copy->dfunc->setSecondaryThresholdDamageDensity(.99);
	
	if(getExtra2dMeshes())
	{
		for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
			copy->addMesh((*getExtra2dMeshes())[i]);
	}
	if(getExtra3dMeshes())
	{
		for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
			copy->addMesh((*getExtra3dMeshes())[i]);
	}
	return copy ; 
}

void WeibullDistributedStiffness::setDamageModel(DamageModel * newmodel)
{
	damageModel = newmodel ;
}

