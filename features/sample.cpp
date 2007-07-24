//
// C++ Implementation: sample
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "sample.h"

using namespace Mu ;

Sample::Sample(Feature * father, double x, double y, double originX, double originY) : Rectangle(x, y, originX, originY), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Rectangle(x+.05, y+.05, originX, originY) ;
	this->boundary2 = new Rectangle(x-.05, y-.05, originX, originY) ;
}

Sample::Sample( double x, double y, double originX, double originY) : Rectangle(x, y, originX, originY), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Rectangle(x+.05, y+.05, originX, originY) ;
	this->boundary2 = new Rectangle(x-.05, y-.05, originX, originY) ;
}

Point * Sample::pointAfter(size_t i)
{
	Point * to_insert = new Point(*(*this->boundingPoints)[i]*0.5 + *(*this->boundingPoints)[(i+1)%this->boundingPoints->size()]*0.5) ;
	std::valarray<Point *> temp(this->boundingPoints->size()+1) ;
	std::copy(&(*this->boundingPoints)[0], &(*this->boundingPoints)[i], &temp[0]) ;
	temp[i+1] = to_insert ;
	std::copy(&(*this->boundingPoints)[i+1], &(*this->boundingPoints)[this->boundingPoints->size()], &temp[i+2]) ;
	this->boundingPoints->resize(temp.size()) ;
	std::copy(&temp[0],&temp[temp.size()] , &(*this->boundingPoints)[0]) ;
	return to_insert ;
}

std::vector<DelaunayTriangle *> Sample::getTriangles( DelaunayTree * dt) 
{
	std::vector<DelaunayTriangle *>  ret ;
	std::vector<DelaunayTriangle *>  temp ;
	if(this->m_f == NULL)
		temp = dt->getTriangles() ;
	else
		temp = dt->conflicts(dynamic_cast<const Rectangle *>(this->boundary)) ;
	
	for(size_t i = 0 ; i < temp.size() ; i++)
	{
		bool inChild = false ;
		for(size_t j = 0 ;  j< this->getChildren()->size() ;  j++)
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

bool Sample::interacts(Feature * f) const
{
	for(Point ** i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary((*i)))
			return true ;
	return false ;
}

