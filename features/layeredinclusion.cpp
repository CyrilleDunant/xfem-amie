//
// C++ Implementation: layeredinclusion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "layeredinclusion.h"


namespace Amie {

std::vector<DelaunayTriangle *> LayeredInclusion::getElements2D( FeatureTree * dt) 
{
	std::vector<DelaunayTriangle *>ret;
	
	std::vector<DelaunayTriangle *>temp = dt->get2DMesh()->getConflictingElements(this->getPrimitive()) ;
	
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
		if(this->in(temp[i]->getCenter()) && inChild == false)
			ret.push_back(temp[i]) ;
	}
	return ret ;
}

LayeredInclusion::LayeredInclusion(Feature *father, std::vector<double> radii,double x,double y) : CompositeFeature(father), LayeredCircle(radii, x, y )
{
	
	for(int i = getRadii().size()-1 ; i > -1 ; i--)
	{
		getComponents().push_back(new VirtualLayer(this, getRadii()[i], getCenter() ) );
	}
}

LayeredInclusion::LayeredInclusion(Feature *father,std::vector<double> radii, const Point center) : CompositeFeature(father), LayeredCircle(radii, center )
{
	
	for(int i = getRadii().size()-1 ; i > -1 ; i--)
	{
		getComponents().push_back(new VirtualLayer(this, getRadii()[i], getCenter() ) );
	}
}

LayeredInclusion::LayeredInclusion(std::vector<double> r,double x,double y) :  CompositeFeature(nullptr), LayeredCircle(r, x, y)
{
	
	for(int i = getRadii().size()-1 ; i > -1 ; i--)
	{
		getComponents().push_back(new VirtualLayer(this, getRadii()[i], getCenter() ) );
	}
}

void LayeredInclusion::print() const
{
	std::cout << "I am a layered inclusion, radius " << getRadius() << ", center (" << getCenter().getX() << "; " << getCenter().getY() << ")"<< std::endl ;
}

LayeredInclusion::LayeredInclusion(std::vector<double> r,Point center) : CompositeFeature(nullptr), LayeredCircle(r, center)
{
	for(int i = getRadii().size()-1 ; i > -1 ; i--)
	{
		getComponents().push_back(new VirtualLayer(this, getRadii()[i], getCenter() ) );
	}
}

LayeredInclusion::LayeredInclusion(double r, Point center) : CompositeFeature(nullptr), LayeredCircle(r, center)
{
	this->isEnrichmentFeature = false ;
	for(int i = getRadii().size()-1 ; i > -1 ; i--)
	{
		getComponents().push_back(new VirtualLayer(this, getRadii()[i], getCenter() ) );
	}
}


void LayeredInclusion::sample(double linearDensity, double surfaceDensityFactor)
{
	this->sampleSurface(linearDensity, surfaceDensityFactor) ;

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

bool LayeredInclusion::interacts(Feature * f, double d) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}


Form * LayeredInclusion::getBehaviour(const Point & p) const
{	
	double pRadius = dist(p, getCenter()) ;
	
	if(layeredBehaviour.size() == 1)
		return behaviour ;
	
	for(size_t i = 0 ; i < getRadii().size() ; i++)
	{
		if (pRadius < getRadii()[i])
			return layeredBehaviour[i] ;
	}
	
	return behaviour ;
}

void LayeredInclusion::setBehaviour(Form * b)
{

	for(size_t i = 0 ; i < layeredBehaviour.size() ; i++)
	{
		if(layeredBehaviour[i] == behaviour)
			layeredBehaviour[i] = nullptr ;
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
			layeredBehaviour[i] = nullptr ;
		else
			delete layeredBehaviour[i] ;
	}

	delete behaviour ;
	behaviour = b[0] ;
	layeredBehaviour.clear() ;
	
	for(size_t i = 0 ; i < std::min(b.size(), getRadii().size()) ; i++)
	{
		layeredBehaviour.push_back(b[i]) ;
	}
}



VirtualLayer::VirtualLayer(LayeredInclusion *father, double r, double x, double y) : VirtualFeature(father), Circle(r, x, y)
{
	source = father ;
}

VirtualLayer::VirtualLayer(LayeredInclusion *father, double r,  Point center) : VirtualFeature(father), Circle(r, center )
{
	source = father ;
}


bool VirtualLayer::interacts(Feature * f, double d) const 	
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}

void VirtualLayer::print() const
{
	std::cout << "I am a virtual layer, radius " << getRadius() << ", center (" << getCenter().getX() << "; " << getCenter().getY() << ")"<< std::endl ;
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

std::vector<DelaunayTriangle *> VirtualLayer::getElements2D( FeatureTree * dt)  { 
	std::vector<DelaunayTriangle *> ret  ;
	
	std::vector<DelaunayTriangle *> temp = dt->get2DMesh()->getConflictingElements(dynamic_cast<Circle *>(this)) ;
	
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
		if(this->in(temp[i]->getCenter()) || inChild )
			ret.push_back(temp[i]) ;
	}
	return ret ;
}

std::vector<DelaunayTetrahedron *> VirtualLayer::getElements3D( FeatureTree* dt)  {
	return std::vector<DelaunayTetrahedron *>(0)  ;
}

Form * VirtualLayer::getBehaviour(const Point & p) const
{
	return source->getBehaviour(p) ;
}

void VirtualLayer::sample(double linearDensity, double surfaceDensityFactor)
{

}

Feature * VirtualLayer::getSource()
{
	return source ;
}

}

