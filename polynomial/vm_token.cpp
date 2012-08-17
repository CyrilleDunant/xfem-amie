//
// C++ Implementation: vm_token
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "vm_token.h"
#include "../elements/elements.h"

using namespace Mu ;


void InHomogeneousProjectionOperation::eval(double * a, double * b, double * c) const
{
	double res = 0;
	Point test(*a, *b) ;
	Point gProj(test) ;
	inGeo->project(&gProj);
	if(inGeo->in(test))
	{
		double totaldist ;
		std::vector<double> weight ;
		std::vector<Point> projs ;
		for(size_t i = 0 ; i < inProjector.size() ; i++)
		{
			Line l(inProjector[i].first(), inProjector[i].vector()) ;
			
			Point proj = l.projection(test) ;

			projs.push_back(proj);
			double d = dist(proj, test) ;
			totaldist += d ;
			weight.push_back(d);
		}
		projs.push_back(gProj);
		double d = dist(gProj, test) ;
		totaldist += d ;
		weight.push_back(d);
		
		double maxn = 0 ;
		double renorm = 0 ;
		for(size_t i = 0 ; i < inProjector.size() ; i++)
		{
			double n = inProjector[i].norm() ;
			if(weight[i] < POINT_TOLERANCE_2D)
			{
				*c = dist(projs[i], inProjector[i].second())/n ;
				return ;
			}
			maxn = std::max(maxn, n) ;
			res += (1./weight[i])*dist(projs[i], inProjector[i].second())/n ;
			renorm += 1./weight[i] ;
		}
		if(weight.back() < POINT_TOLERANCE_2D)
		{
			*c = 1 ;
			return ;
		}
		res += (1./weight.back())*(1.-dist(projs.back(), test)/inGeo->getRadius()) ;
		renorm += 1./weight.back() ;
		
		if(renorm > POINT_TOLERANCE_2D)
			res /= renorm ;
		else
			res = 0 ;
	}
	else
	{
		double totaldist ;
		std::vector<double> weight ;
		std::vector<Point> projs ;
		for(size_t i = 0 ; i < outProjector.size() ; i++)
		{
			Line l(outProjector[i].first(), outProjector[i].vector()) ;
			Point proj = l.projection(test) ;
			projs.push_back(proj);
			double d = dist(proj, test) ;
			totaldist += d ;
			weight.push_back(d);
		}
		projs.push_back(gProj);
		double d = dist(gProj, test) ;
		totaldist += d ;
		weight.push_back(d);
		
		double maxn = 0 ;
		double renorm = 0 ;
		for(size_t i = 0 ; i < outProjector.size() ; i++)
		{
			double n = outProjector[i].norm()  ;
			if(weight[i] < POINT_TOLERANCE_2D)
			{
				*c = dist(projs[i], outProjector[i].second())/n ;
				return ;
			}
			maxn = std::max(maxn, n) ;
			res += (1./weight[i])*dist(projs[i], outProjector[i].second())/n ;
			renorm += 1./weight[i] ;
		}
		if(weight.back() < POINT_TOLERANCE_2D)
		{
			*c = 1 ;
			return ;
		}
		
		res += (1./weight.back())*(1.-dist(projs.back(), test)/inGeo->getRadius()) ;
		renorm += 1./weight.back() ;
		
		if(renorm > POINT_TOLERANCE_2D)
			res /= renorm ;
		else
			res = 0 ;
	}
	
	*c =  res;
}
	

double sign(const double t)
{
	if(t < 0)
		return -1 ;
	if(t > 0)
		return 1 ;
	return 0 ;
} 

double positivity(const double t)
{
	return t >= 0 ;
} 

double negativity(const double t)
{
	return t <= 0 ;
} 

double interpolate(const double a, const double b)
{	
	return a/(b + a) ;
}

