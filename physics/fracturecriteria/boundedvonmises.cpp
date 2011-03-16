
//
// C++ Implementation: vonmises
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "boundedvonmises.h"
namespace Mu {

BoundedVonMises::BoundedVonMises(double thres, double damageThreshold, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, threshold(thres), dmodel(NULL), damageThreshold(damageThreshold)
{

}


BoundedVonMises::~BoundedVonMises()
{
}

double BoundedVonMises::grade(const ElementState &s)
{
	dmodel = s.getParent()->getBehaviour()->getDamageModel() ;
	if(dmodel && dmodel->getState().max() > damageThreshold)
		return -1. ;
		
	double maxStress = s.getMaximumVonMisesStress() ;
	
	if(maxStress > threshold )
	{
		return 1. - std::abs(threshold/maxStress) ;
	}
	else 
	{
		return -1.+ std::abs(maxStress/threshold);
	}
}

FractureCriterion * BoundedVonMises::getCopy() const
{
	return new BoundedVonMises(threshold, damageThreshold) ;
}

Material BoundedVonMises::toMaterial()
{
	Material mat ;
	return mat ;
}



}
