//
// C++ Implementation: sample3d
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2006-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "sample3d.h"

using namespace Amie ;



Sample3D::Sample3D(Feature *father, double x, double y, double z, double originX, double originY, double originZ) : Hexahedron(x, y, z, originX, originY, originZ), Feature(father)
{
	this->isEnrichmentFeature = false ;
}
	
Sample3D::Sample3D(double x, double y, double z,double originX, double originY, double originZ): Hexahedron(x, y, z, originX, originY, originZ), Feature(nullptr)
{
	this->isEnrichmentFeature = false ;
}
	
bool Sample3D::interacts(Feature * f, double d) const {
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}
		
std::vector<DelaunayTriangle *> Sample3D::getElements2D( FeatureTree * dt)
{
	return std::vector<DelaunayTriangle *>(0) ;
}
	
std::vector<DelaunayTetrahedron *> Sample3D::getElements3D( FeatureTree * dt)
{
	std::vector<DelaunayTetrahedron *> ret ;
	std::vector<DelaunayTetrahedron *> temp = dt->get3DMesh()->getConflictingElements(dynamic_cast<const Hexahedron *>(this->getPrimitive())) ;
	
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
		if(this->Hexahedron::in(temp[i]->getCenter()) && inChild == false)
			ret.push_back(temp[i]) ;
	}
	
	return ret ;
}
