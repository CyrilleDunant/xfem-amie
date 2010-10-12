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
#include "physics_base.h"
#include "stiffness.h"
#include "stiffness_and_fracture.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "fracturecriteria/confinedmohrcoulomb.h"
#include "fracturecriteria/mcft.h"
#include "fracturecriteria/confinedmohrcoulombwithstrain.h"
#include "fracturecriteria/maxstrain.h"
#include "fracturecriteria/ruptureenergy.h"
#include "../utilities/random.h"

using namespace Mu ;

WeibullDistributedStiffness::WeibullDistributedStiffness(const Matrix & rig, double down, double up) : LinearForm(rig, true, true, rig.numRows()/3+1), variability(.1), up(up), down(down)
{
	materialRadius = .001;
	neighbourhoodRadius = .004 ;
		
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);

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
	StiffnessAndFracture * ret = new StiffnessAndFracture(
															param*factor, 
															new MCFT(
																					up*factor,
																					down*factor 
																					)
														 ) ;
	ret->criterion->setMaterialCharacteristicRadius(materialRadius);
	ret->criterion->setNeighbourhoodRadius(neighbourhoodRadius);
	ret->dfunc->setMaterialCharacteristicRadius(materialRadius);
	ret->dfunc->setThresholdDamageDensity(.2);
	ret->dfunc->setSecondaryThresholdDamageDensity(.2);
// 	ret->dfunc->setThresholdDamageDensity(.9999);
// 	ret->dfunc->setSecondaryThresholdDamageDensity(.9999);
	return ret ;
}


