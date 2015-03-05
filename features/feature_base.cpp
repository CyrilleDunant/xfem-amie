// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "feature_base.h"
#include "../physics/void_form.h"
#include "../mesher/delaunay.h"
#include "../mesher/delaunay_3d.h"
#include "features.h"

using namespace Amie ;


Feature::Feature(Feature * father) :  m_f(father), behaviour(nullptr)
{
	isEnrichmentFeature = false ;
	isCompositeFeature = false ;
	isVirtualFeature = false ;
	isUpdated = false ;
	behaviourSource = nullptr ;
	
	layer = -1 ;
	if(father != nullptr)
		father->addChild(this) ;

}

std::vector<Point *> Feature::doubleSurfaceSampling()
{
	std::vector<Point *> ret ;
	Amie::PointArray newboundingPoints(getBoundingPoints().size()*2) ;
	for(size_t i = 0 ; i < getBoundingPoints().size()-1 ; i++)
	{
		newboundingPoints[i*2] = &getBoundingPoint(i) ;
		Point * p = new Point(getBoundingPoint(i)*0.5 + getBoundingPoint(i+1)*0.5) ;
		project(p) ;
		p->setId(-1) ;
		newboundingPoints[i*2+1] = p ;
		ret.push_back(newboundingPoints[i*2+1]) ;
	}
	
	newboundingPoints[(getBoundingPoints().size()-1)*2] = &getBoundingPoint(getBoundingPoints().size()-1) ;
	newboundingPoints[getBoundingPoints().size()*2-1] = new Point(getBoundingPoint(getBoundingPoints().size()-1)*0.5 + getBoundingPoint(0)*0.5) ;
	newboundingPoints[getBoundingPoints().size()*2-1]->setId(-1) ;
	project(newboundingPoints[getBoundingPoints().size()*2-1]) ;
	ret.push_back(newboundingPoints[getBoundingPoints().size()*2-1]) ;
	
	dynamic_cast<Geometry *>(this)->setBoundingPoints(newboundingPoints) ;
	
	return ret ;
}


const Feature * Feature::getBehaviourSource() const
{
	return behaviourSource ;
}

Feature * Feature::getBehaviourSource() 
{
	return behaviourSource ;
}

void Feature::setBehaviourSource(Feature * const f)
{
	behaviourSource = f ;
}

void  Feature::addChild( Feature * const feat)
{
	if(getChildren().empty() || std::find(getChildren().begin(), getChildren().end(), feat) == getChildren().end())
	{
		getChildren().push_back(feat) ;
	}
}

void  Feature::removeChild(Feature * const feat)
{
	auto c = std::find(m_c.begin(), m_c.end(), feat) ;
	if(c != m_c.end())
		getChildren().erase(c) ;
}

void Feature::setBehaviour(Form * const f)
{
/*	if(behaviour)
		delete this->behaviour ;*/
	behaviour = f ;
}

Form * Feature::getBehaviour( const Point & p) const 
{
	return this->behaviour ;
}

Form * Feature::getBehaviour() const 
{
	return this->behaviour ;
}

Feature * Feature::getChild(size_t i) const
{
	return m_c[i] ;
}

const std::vector<Feature *> & Feature::getChildren() const
{
	return m_c ;
}

std::vector<Feature *> & Feature::getChildren()
{
	return m_c ;
}

bool Feature::inBoundary(const Point &p, double d) const 
{
	if(p == getCenter())
		return true ;
	Point proj(p) ;
	project(&proj) ;
	return (squareDist3D(proj, p) < d*d) || in(p);
}

bool Feature::onBoundary(const Point &p, double d) const 
{
	if(p == getCenter())
		return false ;
	Point proj(p) ;
	project(&proj) ;
	return (dist(proj, p) < d) ;
}

std::vector<Feature *> Feature::getDescendants() const
{
	std::vector<Feature *> childrenToCheck = getChildren() ;
	std::vector<Feature *> ret ;
	
	while(!childrenToCheck.empty())
	{
		std::vector<Feature *> newChildren ;
		
		for(size_t i = 0 ; i< childrenToCheck.size() ;i++)
		{
			ret.push_back(childrenToCheck[i]) ;
			newChildren.insert(newChildren.end(), childrenToCheck[i]->getChildren().begin(), childrenToCheck[i]->getChildren().end()) ;
		}
		
		childrenToCheck = newChildren ;
	}
	
	return ret ;
}

void  Feature::setFather(Feature * const f)
{
	m_f = f ;
}


CompositeFeature::~CompositeFeature()
{
	for(size_t i = 0 ; i < components.size() ; i++)
		delete components[i] ;
}

std::vector<VirtualFeature *> & CompositeFeature::getComponents()
{
	return components ;
}

const std::vector<VirtualFeature *> & CompositeFeature::getComponents() const
{
	return components ;
}

const Feature * Feature::getFather() const
{
	return m_f ;
}

Feature * Feature::getFather() 
{
	return m_f ;
}



Feature::~Feature() 
{ 
// 	for(size_t i = 0 ; i < this->boundary->size() ; i++)
// 		delete this->boundary->getPoint(i) ;
	
	delete this->behaviour ;
}

std::vector<DelaunayTriangle *> Feature::getBoundingElements2D( FeatureTree * dt) const 
{
	std::vector<DelaunayTriangle *> tri = dt->get2DMesh()->getConflictingElements(dynamic_cast<const Geometry *>(this)) ;
	
	std::vector<DelaunayTriangle *> ret  ;
	
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		if(tri[i]->intersects(dynamic_cast<const Geometry *>(this)))
			ret.push_back(tri[i]) ;
	}
	
	return ret ;

}

std::vector<DelaunayTetrahedron *> Feature::getBoundingElements3D( FeatureTree * dt) const
{
	std::vector<DelaunayTetrahedron *> tri = dt->get3DMesh()->getConflictingElements(dynamic_cast<const Geometry *>(this)) ;
	
	std::vector<DelaunayTetrahedron *> ret  ;
	
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		if(tri[i]->Tetrahedron::intersects(dynamic_cast<const Geometry *>(this)))
			ret.push_back(tri[i]) ;
	}
	
	return ret ;
	
}

