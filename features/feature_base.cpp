#include "feature_base.h"
#include "../physics/void_form.h"
#include "../mesher/delaunay.h"
#include "../mesher/delaunay_3d.h"

using namespace Mu ;


void Feature::setBoundary(Geometry * g)
{
	delete boundary ;
	this->boundary = g ;
}

Geometry * Feature::getBoundary()
{
	return this->boundary ;
}

const Geometry * Feature::getBoundary() const
{
	return this->boundary ;
}

Feature::Feature(Feature * father)
{
	this->isEnrichmentFeature = false ;
	this->isCompositeFeature = false ;
	this->isVirtualFeature = false ;
		
	this->behaviour = new VoidForm() ;

	
	m_f = father ;
	if(father != NULL)
		father->addChild(this) ;
	infRad = DEFAULT_BOUNDARY ;
	this->boundary = NULL ;
	this->boundary2 = NULL ;
}

Feature::Feature()
{
	this->isEnrichmentFeature = false ;
	this->isCompositeFeature = false ;
	this->isVirtualFeature = false ;

	this->behaviour = new VoidForm() ;
	
	m_f = NULL ;
	infRad = DEFAULT_BOUNDARY ;
	this->boundary = NULL ;
	this->boundary2 = NULL ;
}

Feature::Feature(Feature *father, Geometry * b) 
{ 
	boundary = b ; 
	m_f = father ; 
	this->isEnrichmentFeature = false ;
	this->isCompositeFeature = false ;
	this->isVirtualFeature = false ;

	this->behaviour = new VoidForm() ;
	
	if(father != NULL)
		father->addChild(this) ;
	
	infRad = DEFAULT_BOUNDARY ;
	this->boundary2 = NULL ;
}

 bool Feature::inBoundary(const Point & v) const 
{
	bool ret(false) ;
	if(boundary)
		ret = getBoundary()->in(v) ;

	for(size_t i = 0 ;  i < this->m_c.size() ; i++)
		ret = ret || m_c[i]->inBoundary(v) ;
	
	return ret ;
}


std::vector<Point *> Feature::doubleSurfaceSampling()
{
	std::vector<Point *> ret ;
	std::valarray<Point *> newboundingPoints(this->getBoundingPoints().size()*2) ;
	for(size_t i = 0 ; i < this->getBoundingPoints().size()-1 ; i++)
	{
		newboundingPoints[i*2] = &getBoundingPoint(i) ;
		newboundingPoints[i*2+1] = new Point(getBoundingPoint(i)*0.5 + getBoundingPoint(i+1)*0.5) ;
		this->project(newboundingPoints[i*2+1]) ;
		newboundingPoints[i*2+1]->id = -1 ;
		ret.push_back(newboundingPoints[i*2+1]) ;
	}
	
	newboundingPoints[(getBoundingPoints().size()-1)*2] = &getBoundingPoint(getBoundingPoints().size()-1) ;
	newboundingPoints[getBoundingPoints().size()*2-1] = new Point(getBoundingPoint(getBoundingPoints().size()-1)*0.5 + getBoundingPoint(0)*0.5) ;
	newboundingPoints[getBoundingPoints().size()*2-1]->id = -1 ;
	this->project(newboundingPoints[getBoundingPoints().size()*2-1]) ;
	ret.push_back(newboundingPoints[getBoundingPoints().size()*2-1]) ;
	
	dynamic_cast<Geometry *>(this)->setBoundingPoints(newboundingPoints) ;
	
	return ret ;
}

bool Feature::inBoundary(const Point *v) const
{
	bool ret = boundary->in(*v) ;
	for(size_t i = 0 ;  i < this->m_c.size() ; i++)
		ret =  ret || m_c[i]->inBoundary(v) ;
	
	return ret ;
}

bool Feature::inBoundaryLayer(const Point *v) const
{

	if(boundary2)
		return  boundary->in(*v) && !(boundary2->in(*v)) ;
	else
		return boundary->in(*v) ;

}

void  Feature::addChild(Feature *f)
{
	if(std::find(m_c.begin(), m_c.end(), f) == m_c.end())
	{
		m_c.push_back(f) ;
	}
}

void  Feature::removeChild(Feature *f)
{
	std::vector<Feature *>::iterator c = std::find(m_c.begin(), m_c.end(), f) ;
	if(c != m_c.end())
		m_c.erase(c) ;
}

void Feature::setBehaviour(Form * f)
{
	delete this->behaviour ;
	this->behaviour = f ;
}

Form * Feature::getBehaviour( const Point & p)
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

std::vector<Feature *> Feature::getDescendants() const
{
	std::vector<Feature *> childrenToCheck = m_c ;
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

void  Feature::setFather(Feature *f)
{
	m_f = f ;
}

 void Feature::setInfluenceRadius(double r)
{
	infRad = r ;
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

Feature * Feature::getFather() const
{
	return this->m_f ;
}


Feature::~Feature() 
{ 
// 	for(size_t i = 0 ; i < this->boundary->size() ; i++)
// 		delete this->boundary->getPoint(i) ;
	
	delete this->boundary ;
	delete this->boundary2 ;
	
	delete this->behaviour ;
}

std::vector<DelaunayTriangle *> Feature::getBoundingElements( Mesh<DelaunayTriangle> * dt)
{
	std::vector<DelaunayTriangle *> tri = dt->getConflictingElements(dynamic_cast<Geometry *>(this)) ;
	
	std::vector<DelaunayTriangle *> ret  ;
	
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		if(tri[i]->intersects(dynamic_cast<Geometry *>(this)))
			ret.push_back(tri[i]) ;
	}
	
	return ret ;

}

std::vector<DelaunayTetrahedron *> Feature::getBoundingElements( Mesh<DelaunayTetrahedron> * dt)
{
	std::vector<DelaunayTetrahedron *> tri = dt->getConflictingElements(dynamic_cast<Geometry *>(this)) ;
	
	std::vector<DelaunayTetrahedron *> ret  ;
	
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		if(tri[i]->Tetrahedron::intersects(dynamic_cast<Geometry *>(this)))
			ret.push_back(tri[i]) ;
	}
	
	return ret ;
	
}

