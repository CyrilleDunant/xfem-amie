//
// C++ Implementation: layeredinclusion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "layeredinclusion.h"


using namespace Mu ;

std::vector<DelaunayTriangle *> LayeredInclusion::getTriangles( DelaunayTree * dt) 
{
	std::vector<DelaunayTriangle *>ret;
	
	std::vector<DelaunayTriangle *>temp = dt->conflicts(this->boundary) ;
	
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
		if(this->boundary->in(temp[i]->getCenter()) && inChild == false)
			ret.push_back(temp[i]) ;
	}
	return ret ;
}

LayeredInclusion::LayeredInclusion(Feature *father, std::vector<double> radii,double x,double y) : LayeredCircle(radii, x, y ), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Circle(getRadius()*1.1, x, y) ;
	this->boundary2 = new Circle(getRadius()*.9, x, y) ;
}

LayeredInclusion::LayeredInclusion(Feature *father,std::vector<double> radii, const Point center) : LayeredCircle(radii, center ), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Circle(getRadius()*1.1, center) ;
	this->boundary2 = new Circle(getRadius()*.9, center) ;
}

LayeredInclusion::LayeredInclusion(std::vector<double> r,double x,double y) : LayeredCircle(r, x, y), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Circle(getRadius()*1.1, x, y) ;
	this->boundary2 = new Circle(getRadius()*.9, x, y) ;
}

LayeredInclusion::LayeredInclusion(std::vector<double> r,Point center) : LayeredCircle(r, center), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Circle(getRadius()*1.1, center) ;
	this->boundary2 = new Circle(getRadius()*.9, center) ;
}

LayeredInclusion::LayeredInclusion(double r, Point center) : LayeredCircle(r, center),  Feature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Circle(getRadius()*1.1, center) ;
	this->boundary2 = new Circle(getRadius()*.9, center) ;
}

Point * LayeredInclusion::pointAfter(size_t i)
{
	double theta_i = atan2(boundingPoints[i]->y, boundingPoints[i]->x) ;
	double theta_ip = atan2(boundingPoints[(i+1)%this->boundingPoints.size()]->y, boundingPoints[(i+1)%this->boundingPoints.size()]->x) ;
	double theta = 0.5*theta_i + 0.5*theta_ip ;
	
	Point * to_insert = new Point(cos(theta)*this->getRadius()+ getCenter().x, sin(theta)*getRadius()+ getCenter().y) ;
	std::valarray<Point *> temp(this->boundingPoints.size()+1) ;
	std::copy(&boundingPoints[0], &boundingPoints[i], &temp[0]) ;
	temp[i+1] = to_insert ;
	std::copy(&boundingPoints[i+1], &boundingPoints[this->boundingPoints.size()], &temp[i+2]) ;
	this->boundingPoints.resize(temp.size()) ;
	std::copy(&temp[0],&temp[temp.size()] , &boundingPoints[0]) ;
	return to_insert ;
}

void LayeredInclusion::sample(size_t n)
{
	delete this->boundary ;
	double numberOfRings =round((double)n/(2. * M_PI )) ;
	double r = getRadius()*(1.+ .6/(numberOfRings+1)) ;
	this->boundary = new Circle(r, getCenter()) ;
	this->sampleSurface(n) ;
}

std::vector<Geometry *> LayeredInclusion::getRefinementZones(size_t level) const
{
	std::vector<Geometry *> ret ;
	if(level > 0)
		ret.push_back(new Circle(getRadius()*2., getCenter())) ;
	if(level > 1)
		ret.push_back(new Circle(getRadius() * 1.5, getCenter())) ;
	if(level > 2)
		ret.push_back(new Circle(getRadius() * 1.1, getCenter())) ;
	return ret ;
}

bool LayeredInclusion::interacts(Feature * f) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary((*i)))
			return true ;
	return false ;
}

// Form * LayeredInclusion::getBehaviour(const Point & p)
// {
// 	//special case for Layered Parent ;
// 	if(m_f && m_f->getGeometryType() == LAYERED_CIRCLE && m_f->in(p))
// 	{
// 		std::vector<double> parentRadii = dynamic_cast<LayeredCircle *>(m_f)->getRadii() ;
// 		std::vector<double> childRadii = getRadii() ;
// 		Point parentCenter = m_f->getCenter() ;
// 		Point childCenter = getCenter() ;
// 		
// 		std::vector<double> rlist ;
// 		std::vector<Point> clist ;
// 		std::vector<bool> inchild ;
// 
// 		for(int i = std::min(parentRadii.size(), childRadii.size())-1 ; i >= 0 ; i--)
// 		{
// 			rlist.push_back(parentRadii[i]) ;
// 			inchild.push_back(false) ;
// 			clist.push_back(parentCenter) ;
// 			rlist.push_back(childRadii[i]) ;
// 			inchild.push_back(true) ;
// 			clist.push_back(childCenter) ;
// 		}
// 		
// 		int index = 0 ;
// 		
// 		for(int i = 1 ; i <rlist.size()  ; i++)
// 		{
// 			index = i ;
// 			if(dist(p, clist[i]) > rlist[i])
// 				break ;
// 		}
// 		
// 		if(!inchild[index])
// 			return dynamic_cast<LayeredInclusion *>(m_f)->getBehaviour(p, true) ;
// 		
// 	}
// 
// 	getBehaviour(p, true) ;
// }

Form * LayeredInclusion::getBehaviour(const Point & p)
{	
	double pRadius = dist(p, getCenter()) ;
	
	if(layeredBehaviour.size() == 1)
		return behaviour ;
	
	for(size_t i = 0 ; i < radiuses.size() ; i++)
	{
		if (pRadius < radiuses[i])
			return layeredBehaviour[i] ;
	}
	
	return behaviour ;
}

void LayeredInclusion::setBehaviour(Form * b)
{

	for(size_t i = 0 ; i < layeredBehaviour.size() ; i++)
	{
		if(layeredBehaviour[i] == behaviour)
			layeredBehaviour[i] = NULL ;
		else
			delete layeredBehaviour[i] ;
	}

	layeredBehaviour.clear() ;

	delete behaviour ;
	behaviour = b ;
	layeredBehaviour.push_back(b) ;

}

void LayeredInclusion::setBehaviours(std::vector<Form *> b)
{
	for(size_t i = 0 ; i < layeredBehaviour.size() ; i++)
	{
		if(layeredBehaviour[i] == behaviour)
			layeredBehaviour[i] = NULL ;
		else
			delete layeredBehaviour[i] ;
	}

	delete behaviour ;
	behaviour = b[0] ;
	layeredBehaviour.clear() ;
	
	for(size_t i = 0 ; i < std::min(b.size(), radiuses.size()) ; i++)
	{
		layeredBehaviour.push_back(b[i]) ;
	}
}



VirtualLayer::VirtualLayer(LayeredInclusion *father, double r, double x, double y) : Circle(r, x, y), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->boundary = NULL ; //new Sphere(r*1.07, x, y,z) ;
	this->boundary2 =  NULL ; //new Sphere(r*0.93, x, y,z) ;
}

VirtualLayer::VirtualLayer(LayeredInclusion *father, double r,  Point center) : Circle(r, center ), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->boundary =  NULL ; //new Sphere(r*1.07, center) ;
	this->boundary2 =  NULL ; //new Sphere(r*0.93, center) ;
}


bool VirtualLayer::interacts(Feature * f) const 	
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary((*i)))
			return true ;
	return false ;
}

bool VirtualLayer::inBoundary(const Point *) const
{
	return false ;
}

bool VirtualLayer::inBoundary(const Point &) const
{
	return false ;
}

std::vector<Geometry *> VirtualLayer::getRefinementZones(size_t level) const 
{
	std::vector<Geometry *> ret ;
	if(level > 0)
		ret.push_back(new Circle(sqrt(radius)*(1.2), Circle::getCenter())) ;
	if(level > 1)
		ret.push_back(new Circle(sqrt(radius)*(1.15), Circle::getCenter())) ;
	if(level > 2)
		ret.push_back(new Circle(sqrt(radius)*(1.08), Circle::getCenter())) ;
	return ret ;
}

std::vector<DelaunayTriangle *> VirtualLayer::getTriangles( DelaunayTree * dt)  { 
	std::vector<DelaunayTriangle *> ret  ;
	
	std::vector<DelaunayTriangle *> temp = dt->conflicts(dynamic_cast<Circle *>(this)) ;
	
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
		if(this->in(temp[i]->getCenter()) || inChild )
			ret.push_back(temp[i]) ;
	}
	return ret ;
}

std::vector<DelaunayTetrahedron *> VirtualLayer::getTetrahedrons( DelaunayTree_3D * dt)  {
	return std::vector<DelaunayTetrahedron *>(0)  ;
}


Point * VirtualLayer::pointAfter(size_t i) { 
	return NULL ; 
}

Form * VirtualLayer::getBehaviour(const Point & p)
{
	return m_f->getBehaviour(p) ;
}

void VirtualLayer::sample(size_t n)
{
}

