//
// C++ Implementation: samplingcriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
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
	Point v0 = Point(t->second->x - t->first->x, t->second->y - t->first->y) ;
	Point v1 = Point(t->third->x - t->second->x, t->third->y - t->second->y) ;
	Point v2 = Point(t->first->x - t->third->x, t->first->y - t->third->y) ;
	std::vector<double> ang;
	ang.push_back(
	               asin(
	                     ( v0.x*v1.y - v0.y*v1.x )/( v0.norm()*v1.norm() )
	                   )
	             ) ;
	ang.push_back(
	               asin(
	                     ( v1.x*v2.y - v1.y*v2.x )/( v1.norm()*v2.norm() )
	                   )
	             ) ;
	ang.push_back(
	               asin(
	                     ( v2.x*v0.y - v2.y*v0.x )/( v2.norm()*v0.norm() )
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

	Segment s0(*t->first, *t->second) ;
	Segment s1(*t->second, *t->third) ;
	Segment s2(*t->third, *t->first) ;
	std::map<double, Segment> sides;
	sides[s0.vector().norm()] = s0 ;
	sides[s1.vector().norm()] = s1 ;
	sides[s2.vector().norm()] = s2 ;
	
	std::vector<Point> ret ;
	
	ret.push_back(sides.rbegin()->second.midPoint()) ;
	ret.push_back((++sides.begin())->second.midPoint()) ;
	
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
	Point v0 = Point(t->second->x - t->first->x, t->second->y - t->first->y) ;
	Point v1 = Point(t->third->x - t->second->x, t->third->y - t->second->y) ;
	Point v2 = Point(t->first->x - t->third->x, t->first->y - t->third->y) ;
	
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
