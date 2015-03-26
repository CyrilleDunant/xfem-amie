//
// C++ Implementation: ruptureenergy
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ruptureenergy.h"

namespace Amie {

RuptureEnergy::RuptureEnergy(double e, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, energy(e)
{
}


RuptureEnergy::~RuptureEnergy()
{
}

double RuptureEnergy::grade(ElementState &s)
{
	Vector pstress(0., s.getParent()->getGaussPoints().gaussPoints.size()*(3+3*(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL))) ;
	Vector pstrain(0., pstress.size()) ;
	s.getField( REAL_STRESS_FIELD, STRAIN_FIELD, s.getParent()->getGaussPoints().gaussPoints, pstress, pstrain, true) ;
	
	pstrain *= pstress ;
	
	Vector E(pstrain.size()/3) ;
	
	for(size_t j = 0 ; j < E.size() ; ++j)
	{
		E[j] = pstrain[j*3] + pstrain[j*3+1] + pstrain[j*3+2] ;
	}
	double enr = VirtualMachine().ieval(E, s.getParent())/s.getParent()->area() ;
	if( enr > energy )
	{
		return 1. - std::abs(energy/enr) ;
	}
	else 
	{
		return -1.+ std::abs(enr/energy);
	}
}

FractureCriterion * RuptureEnergy::getCopy() const
{
	return new RuptureEnergy(*this) ;
}

}
