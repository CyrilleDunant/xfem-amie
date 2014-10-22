//
// C++ Implementation: sample
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2006-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "sample.h"
#include "../physics/void_form.h"

using namespace Amie ;

Sample::Sample(Feature * father, double x, double y, double originX, double originY) : Rectangle(x, y, originX, originY), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = nullptr ;
}

Sample::Sample( double x, double y, double originX, double originY) : Rectangle(x, y, originX, originY), Feature(nullptr)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = nullptr ;
}

std::vector<DelaunayTriangle *> Sample::getElements2D( FeatureTree * dt) 
{
	std::vector<DelaunayTriangle *>  ret ;
	std::vector<DelaunayTriangle *>  temp = dt->get2DMesh()->getConflictingElements(this->getPrimitive()) ;
	
	for(size_t i = 0 ; i < temp.size() ; i++)
	{
		bool inChild = false ;
		for(size_t j = 0 ;  j< this->getChildren().size() ;  j++)
		{
			if(this->getChild(j)->in(temp[i]->getCenter()))
			{
				inChild = true ; 
				break ;
			}
		}
		if(this->Rectangle::in(temp[i]->getCenter()) && inChild == false)
			ret.push_back(temp[i]) ;
	}
	
	return ret ;
}

bool Sample::interacts(Feature * f, double d) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}

