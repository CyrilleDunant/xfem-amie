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

#include "loftedPolygonSample3d.h"
#include "../physics/void_form.h"

using namespace Amie ;

LoftedPolygonalSample3D::LoftedPolygonalSample3D(Feature *father, const std::valarray<Point *> & points, const std::vector<Point> & interpolationPoints) : LoftedPolygonPrism(points,  interpolationPoints), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = nullptr ;
}


std::vector<DelaunayTetrahedron *> LoftedPolygonalSample3D::getElements3D( FeatureTree * dt) 
{
	std::vector<DelaunayTetrahedron *>  ret ;
	std::vector<DelaunayTetrahedron *>  temp = dt->get3DMesh()->getConflictingElements(getPrimitive()) ;
	
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
		if(this->LoftedPolygonalSample3D::in(temp[i]->getCenter()) && inChild == false)
			ret.push_back(temp[i]) ;
	}
	
	return ret ;
}

bool LoftedPolygonalSample3D::interacts(Feature * f, double d) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}

