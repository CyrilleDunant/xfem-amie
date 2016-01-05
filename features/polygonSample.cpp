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

#include "polygonSample.h"
#include "../physics/void_form.h"

namespace Amie {

PolygonalSample::PolygonalSample(Feature *father, const std::valarray<Point *> & points) : Polygon(points), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = nullptr ;
}


std::vector<DelaunayTriangle *> PolygonalSample::getElements2D( FeatureTree * dt) 
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
		if(this->PolygonalSample::in(temp[i]->getCenter()) && inChild == false)
			ret.push_back(temp[i]) ;
	}
	
	return ret ;
}

bool PolygonalSample::interacts(Feature * f, double d) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}

}