//
// C++ Implementation: samplingcriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
//

#include "samplingcriterion.h"

using namespace Mu ;

SamplingCriterion::SamplingCriterion()
{
	toAdd = true ;
}


bool SamplingCriterion::add() const
{
	return toAdd ;
}

SamplingCriterion::~SamplingCriterion()
{
}

bool SamplingCriterion::meetsCriterion(const DelaunayTriangle * t)  
{
	return true ;
}

MinimumAngle::MinimumAngle(double a) 
{ 
	angle = a ;
}
	
MinimumAngle::~MinimumAngle() 
{
	
}
	
bool MinimumAngle::meetsCriterion(const DelaunayTriangle * t)  
{
	Point v0 = Point(t->second->getX() - t->first->getX(), t->second->getY() - t->first->getY()) ;
	Point v1 = Point(t->third->getX() - t->second->getX(), t->third->getY() - t->second->getY()) ;
	Point v2 = Point(t->first->getX() - t->third->getX(), t->first->getY() - t->third->getY()) ;
	std::vector<double> ang;
	ang.push_back(
	               asin(
	                     ( v0.getX()*v1.getY() - v0.getY()*v1.getX() )/( v0.norm()*v1.norm() )
	                   )
	             ) ;
	ang.push_back(
	               asin(
	                     ( v1.getX()*v2.getY() - v1.getY()*v2.getX() )/( v1.norm()*v2.norm() )
	                   )
	             ) ;
	ang.push_back(
	               asin(
	                     ( v2.getX()*v0.getY() - v2.getY()*v0.getX() )/( v2.norm()*v0.norm() )
	                   )
	             ) ;
	
	std::sort(ang.begin(), ang.end()) ;

	if(ang[0]+ang[1]+ang[2] > 1.000001*M_PI)
	{
		std::cout << "what a wierd triangle..." ; t->print() ;
		std::cout << "angle sum is " << ang[0]+ang[1]+ang[2] << std::endl ;
		assert(false) ;
	}

	return ang[0] > angle ;
}

void MinimumAngle::reset() 
{
}

std::vector<Point> MinimumAngle::suggest(const DelaunayTriangle * t) const
{

	std::vector<Point> ret ;

	double d0 = dist(t->commonEdge(t->getNeighbour(0)).first, t->commonEdge(t->getNeighbour(0)).second) ;
	double d1 = dist(t->commonEdge(t->getNeighbour(1)).first, t->commonEdge(t->getNeighbour(1)).second) ;
	double d2 = dist(t->commonEdge(t->getNeighbour(2)).first, t->commonEdge(t->getNeighbour(2)).second) ;
	if(d0 < d1 && d0 < d2)
	{
		if(t->getNeighbour(1)->isTriangle)
			ret.push_back(dynamic_cast<DelaunayTriangle *>(t->getNeighbour(1))->getCenter()) ;
		else
			ret.push_back(*t->commonEdge(t->getNeighbour(1)).first*.5+ *t->commonEdge(t->getNeighbour(1)).second*.5) ;
		if(t->getNeighbour(2)->isTriangle)
			ret.push_back(dynamic_cast<DelaunayTriangle *>(t->getNeighbour(2))->getCenter()) ;
		else
			ret.push_back(*t->commonEdge(t->getNeighbour(2)).first*.5+ *t->commonEdge(t->getNeighbour(2)).second*.5) ;
	}

	if(d2 < d1 && d2 < d0)
	{
		if(t->getNeighbour(1)->isTriangle)
			ret.push_back(dynamic_cast<DelaunayTriangle *>(t->getNeighbour(1))->getCenter()) ;
		else
			ret.push_back(*t->commonEdge(t->getNeighbour(1)).first*.5+ *t->commonEdge(t->getNeighbour(1)).second*.5) ;
		if(t->getNeighbour(0)->isTriangle)
			ret.push_back(dynamic_cast<DelaunayTriangle *>(t->getNeighbour(0))->getCenter()) ;
		else
			ret.push_back(*t->commonEdge(t->getNeighbour(0)).first*.5+ *t->commonEdge(t->getNeighbour(0)).second*.5) ;
	}
	
	if(d1 < d2 && d1 < d0)
	{
		if(t->getNeighbour(0)->isTriangle)
			ret.push_back(dynamic_cast<DelaunayTriangle *>(t->getNeighbour(0))->getCenter()) ;
		else
			ret.push_back(*t->commonEdge(t->getNeighbour(0)).first*.5+ *t->commonEdge(t->getNeighbour(0)).second*.5) ;
		if(t->getNeighbour(2)->isTriangle)
			ret.push_back(dynamic_cast<DelaunayTriangle *>(t->getNeighbour(2))->getCenter()) ;
		else
			ret.push_back(*t->commonEdge(t->getNeighbour(2)).first*.5+ *t->commonEdge(t->getNeighbour(2)).second*.5) ;
	}
	
	return ret ;
// 	std::vector<Point *> points ;
// 
// 	for(size_t i = 0 ; i < t->neighbourhood.size() ; i++)
// 	{
// 		points.push_back(t->neighbourhood[i]->first) ;
// 		points.push_back(t->neighbourhood[i]->second) ;
// 		points.push_back(t->neighbourhood[i]->third) ;
// 	}
// 	
// 	std::sort(points.begin(), points.end()) ;
// 	std::vector<Point *>::iterator e = std::unique(points.begin(), points.end()) ;
// 	points.erase(e, points.end()) ;
// 	
// 	Point center;
// 	
// 	double xmin = center.getX();
// 	double xmax = center.getX();
// 	double ymin = center.getY();
// 	double ymax = center.getY();
// 	for(size_t i = 0 ; i < points.size() ; i++)
// 	{
// 		center+= *points[i]/points.size() ;
// 		xmin = std::min(xmin, points[i]->getX()) ;
// 		xmax = std::max(xmax, points[i]->getX()) ;
// 		ymin = std::min(ymin, points[i]->getY()) ;
// 		ymax = std::max(ymax, points[i]->getY()) ;
// 	}
// 	
// 	double radiusx = xmax-xmin ;
// 	double radiusy = ymax-ymin ;
// 
// 	std::vector<Point> fixPoints ;
// 	for(size_t i = 0 ; i < points.size() ; i++)
// 	{
// 		
// 		fixPoints.push_back(*points[i]+Point(-radiusx, radiusy)) ;
// 		fixPoints.push_back(*points[i]+Point(0, radiusy)) ;
// 		fixPoints.push_back(*points[i]+Point(radiusx, radiusy)) ;
// 		fixPoints.push_back(*points[i]+Point(-radiusx, 0)) ;
// 		fixPoints.push_back(*points[i]) ;
// 		fixPoints.push_back(*points[i]+Point(radiusx, 0)) ;
// 		fixPoints.push_back(*points[i]+Point(-radiusx, -radiusy)) ;
// 		fixPoints.push_back(*points[i]+Point(0, -radiusy)) ;
// 		fixPoints.push_back(*points[i]+Point(radiusx, -radiusy)) ;
// 	}
// 	
// 	std::vector<Point> newPoints ;
// 	newPoints.push_back(t->getCenter()) ;
// 	for(size_t i = 0 ; i < t->neighbourhood.size() ; i++)
// 		newPoints.push_back(t->neighbourhood[i]->getCenter()) ;
// 	
// 	
// 	std::vector<Point> speed(newPoints.size()) ;
// 	for(size_t i = 0 ; i < 2000 ; i++)
// 	{
// 		std::vector<Point> forces(newPoints.size()) ;
// 		
// 		for(size_t j = 0 ; j < fixPoints.size() ; j++)
// 		{
// 			for(size_t k = 0 ; k < newPoints.size() ; k++)
// 			{
// 				Point dir = fixPoints[j]-newPoints[k] ;
// 				dir /= -dir.norm() ;
// 				forces[k]+=dir / dist(newPoints[k],fixPoints[j]) ;
// 			}
// 		}
// 		
// 		std::vector<Point> mirrorPoints = newPoints ;
// 		for(size_t j = 0 ; j < newPoints.size() ; j++)
// 		{
// 			mirrorPoints.push_back(newPoints[j]+Point(-radiusx, radiusy)) ;
// 			mirrorPoints.push_back(newPoints[j]+Point(0, radiusy)) ;
// 			mirrorPoints.push_back(newPoints[j]+Point(radiusx, radiusy)) ;
// 			mirrorPoints.push_back(newPoints[j]+Point(-radiusx, 0)) ;
// 			mirrorPoints.push_back(newPoints[j]+Point(radiusx, 0)) ;
// 			mirrorPoints.push_back(newPoints[j]+Point(-radiusx, -radiusy)) ;
// 			mirrorPoints.push_back(newPoints[j]+Point(0, -radiusy)) ;
// 			mirrorPoints.push_back(newPoints[j]+Point(radiusx, -radiusy)) ;
// 		}
// 		
// 		for(size_t j = 0 ; j < newPoints.size() ; j++)
// 		{
// 			for(size_t k = j+1 ; k < newPoints.size() ; k++)
// 			{
// 				Point dir = mirrorPoints[k]-newPoints[j] ;
// 				dir /= -dir.norm() ;
// 				forces[j]+=dir / dist(mirrorPoints[k],newPoints[j]) ;
// 			}
// 		}
// 		
// 		for(size_t j = 0 ; j < newPoints.size() ; j++)
// 		{
// 			for(size_t k = newPoints.size() ; k < mirrorPoints.size() ; k++)
// 			{
// 				Point dir = mirrorPoints[k]-newPoints[j] ;
// 				dir /= -dir.norm() ;
// 				forces[j]+=dir / dist(mirrorPoints[k],newPoints[j]) ;
// 			}
// 		}
// 		
// 		for(size_t j = 0 ; j < newPoints.size() ; j++)
// 		{
// 			if(speed[j].norm() > .001)
// 				speed[j] -= speed[j]/speed[j].norm()*.1 ;
// 			speed[j] += forces[j]*(radiusx*radiusx+radiusy*radiusy)*.001 ;
// 		}
// 		for(size_t j = 0 ; j < newPoints.size() ; j++)
// 		{
// 			newPoints[j]+=speed[j] ;
// 			
// 			while(newPoints[j].getX() >= xmax )
// 			{
// 				newPoints[j].getX() -= radiusx ;
// 			}
// 			while(newPoints[j].getX() <= xmin)
// 			{
// 				newPoints[j].getX() += radiusx ;
// 			}
// 			while(newPoints[j].getY() >= ymax)
// 			{
// 				newPoints[j].getY() -= radiusy ;
// 			}
// 			while(newPoints[j].getY() <= ymin)
// 			{
// 				newPoints[j].getY() += radiusy ;
// 			}
// 		}
// 		
// 	}
// 	
// 	for(size_t i = 0 ; i < newPoints.size() ; i++)
// 	{
// 		bool in =t->in(newPoints[i]);
// 		
// 		for(size_t j = 0 ; j < t->neighbourhood.size() ; j++)
// 		{
// 			if(in)
// 			{
// 				ret.push_back(newPoints[i]) ;
// 				break ;
// 			}
// 			
// 			if(t->neighbourhood[j]->in(newPoints[i]))
// 			{
// 				ret.push_back(newPoints[i]) ;
// 				break ;
// 			}
// 		}
// 	}
	
	return ret ;
}

MaximumLength::MaximumLength(double l)
{
	length = l ;
}
	
MaximumLength::~MaximumLength()
{
}
	
bool MaximumLength::meetsCriterion(const DelaunayTriangle * t) 
{
	Point v0 = Point(t->second->getX() - t->first->getX(), t->second->getY() - t->first->getY()) ;
	Point v1 = Point(t->third->getX() - t->second->getX(), t->third->getY() - t->second->getY()) ;
	Point v2 = Point(t->first->getX() - t->third->getX(), t->first->getY() - t->third->getY()) ;
	
	std::vector<double> len;
	len.push_back(v0.norm()) ;
	len.push_back(v1.norm()) ;
	len.push_back(v2.norm()) ;
	
	std::sort(len.begin(), len.end()) ;
	
	return len[0] < length ;
}

void MaximumLength::reset() 
{
}

Counter::Counter(size_t nit)
{
	this->c = 0 ;
	this->s = nit ;
	toAdd = false ;
}
	
Counter::~Counter()
{
}
	
bool Counter::meetsCriterion(const DelaunayTriangle * t) 
{
	c++ ;
	
	return c < s ;
}

void Counter::reset()
{
	c = 0 ;
}
