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

	std::vector<Point> ret ;
	
	std::vector<Point *> points ;
	double maxX = std::max(std::max(t->first->x, t->second->x), t->third->x) ;
	double minX = std::min(std::min(t->first->x, t->second->x), t->third->x) ;
	double maxY = std::max(std::max(t->first->y, t->second->y), t->third->y) ;
	double minY = std::min(std::min(t->first->y, t->second->y), t->third->y) ;
	double count = 1 ;
	double averageLength = (dist(t->first,t->second) + dist(t->first,t->third)+ dist(t->second,t->third))/3.;
	for(size_t i = 0 ; i < t->neighbourhood.size() ; i++)
	{
		points.push_back(t->neighbourhood[i]->first) ;
		points.push_back(t->neighbourhood[i]->second) ;
		points.push_back(t->neighbourhood[i]->third) ;
		averageLength +=  (dist(t->neighbourhood[i]->first,t->neighbourhood[i]->second) 
		                   + dist(t->neighbourhood[i]->first,t->neighbourhood[i]->third)
		                   + dist(t->neighbourhood[i]->second,t->neighbourhood[i]->third))/3.;
		count++ ;
		maxX = std::max(
		                 std::max(
		                           std::max(
		                                     t->neighbourhood[i]->first->x, 
		                                     t->neighbourhood[i]->second->x
		                                   ), 
		                           t->neighbourhood[i]->third->x
		                         ),
		                 maxX
		               ) ;
		minX =  std::min(
		                         std::min(
		                                   std::min(
		                                             t->neighbourhood[i]->first->x, 
		                                             t->neighbourhood[i]->second->x
		                                           ), 
		                                   t->neighbourhood[i]->third->x
		                                 ),
		                         minX
		                       ) ;
		maxY = std::max(
		                        std::max(
		                                  std::max(
		                                            t->neighbourhood[i]->first->y, 
		                                            t->neighbourhood[i]->second->y
		                                          ), 
		                                  t->neighbourhood[i]->third->y
		                                ),
		                        maxY
		                      ) ;
		minY =  std::min(
		                         std::min(
		                                   std::min(
		                                             t->neighbourhood[i]->first->y, 
		                                             t->neighbourhood[i]->second->y
		                                           ), 
		                                   t->neighbourhood[i]->third->y
		                                 ),
		                         minY
		                       ) ;
	}
	
	std::sort(points.begin(), points.end()) ;
	std::vector<Point *>::iterator e = std::unique(points.begin(), points.end()) ;
	points.erase(e, points.end()) ;
	averageLength /= count ;
	int pointsAlongX = ceil((maxX-minX)/averageLength);
	int pointsAlongY = ceil((maxY-minY)/averageLength);
	std::cout << "grid is " << pointsAlongX << "×" << pointsAlongY << std::endl ;
	double d = .5*averageLength ;
	for(int i = 0 ; i < pointsAlongX+2 ; i++)
	{
		for(int j = 0 ; j < pointsAlongY+2 ; j++)
		{
			Point potential(minX+(double)i/(pointsAlongX+1)*(maxX-minX), minY+(double)j/(pointsAlongY+1)*(maxY-minY)) ;
			
			bool inNeighbourhood = false ;
			
			for(size_t k = 0 ; k < t->neighbourhood.size() ; k++)
			{
				if(t->neighbourhood[k]->in(potential))
				{
					inNeighbourhood = true ;
					break ;
				}
			}
			if(t->in(potential))
				inNeighbourhood = true ;
			
			if(inNeighbourhood)
			{
				bool alone = true ;
				for(size_t k = 0 ; k < points.size() ;k++)
				{
					if(dist(points[k], &potential) < d)
					{
						alone = false ;
						break ;
					}
				}
				
				if(alone)
					ret.push_back(potential) ;
			}
		}
	}
// 	ret.push_back(sides.rbegin()->second.midPoint()) ;
// 	ret.push_back(sides.rbegin()->second.midPoint() + ((++sides.begin())->second.midPoint()-sides.rbegin()->second.midPoint()) * 2. ) ;
	
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
