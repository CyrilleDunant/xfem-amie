//
// C++ Implementation: ruptureenergy
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ruptureenergy.h"
#include "../mesher/delaunay.h"

namespace Mu {

RuptureEnergy::RuptureEnergy(double e)
	: energy(e)
{
}


RuptureEnergy::~RuptureEnergy()
{
}

double RuptureEnergy::grade(const ElementState &s) const
{
	Vector pstress = s.getStress(s.getParent()->getGaussPoints().gaussPoints) ;
	Vector pstrain = s.getStrain(s.getParent()->getGaussPoints().gaussPoints) ;
	
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
		return 0 ;
	}
}

FractureCriterion * RuptureEnergy::getCopy() const
{
	return new RuptureEnergy(*this) ;
}

}
