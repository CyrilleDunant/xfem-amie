
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
#include "nonlocalvonmises.h"
#include "../../mesher/delaunay.h"
namespace Mu {

NonLocalVonMises::NonLocalVonMises(double thresh, double radius) : threshold(thresh), radius(radius)
{
}


NonLocalVonMises::~NonLocalVonMises()
{
}

double NonLocalVonMises::grade(const ElementState &s)
{
	if(cache.empty())
	{
		Circle epsilon(radius,s.getParent()->getCenter()) ;
		std::vector<DelaunayTriangle *> tempcache = s.getParent()->get2DMesh()->getConflictingElements(&epsilon);
		for(size_t i = 0 ; i < tempcache.size() ; i++)
		{
			if(tempcache[i]->getBehaviour()->getFractureCriterion())
				cache.push_back(tempcache[i]);
		}
	}
	double maxStress = 0 ;
	double area = 0 ;
	for(int i = 0 ; 0 < cache.size() ; i++)
	{
		double tmparea = cache[i]->area() ;
		maxStress += cache[i]->getState().getMaximumVonMisesStress()*tmparea ;
		area += tmparea ;
	}

	maxStress /= area ;
	
	if(maxStress > threshold )
	{
		return 1. - std::abs(threshold/maxStress) ;
	}
	else 
	{
		return -1.+ std::abs(maxStress/threshold);
	}
}

FractureCriterion * NonLocalVonMises::getCopy() const
{
	return new NonLocalVonMises(threshold, radius) ;
}

Material NonLocalVonMises::toMaterial()
{
	Material mat ;
	mat(TAG_MAX_STRESS,threshold) ;
	return mat ;
}

}
