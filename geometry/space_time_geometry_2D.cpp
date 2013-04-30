#include "space_time_geometry_2D.h"

using namespace Mu ;

TimeDependentCircle::TimeDependentCircle(Function r, const Point & c): Circle(0., c.x, c.y), radius_t(r)
{
	gType = TIME_DEPENDENT_CIRCLE ;
	center.t = c.t ;
}


TimeDependentCircle::TimeDependentCircle(Function r, Point * c): Circle(0., c->x, c->y), radius_t(r)
{
	gType = TIME_DEPENDENT_CIRCLE ;
	center.t = c->t ;
}

TimeDependentCircle::TimeDependentCircle(double r0, double rate, const Point & c): Circle(r0, c.x, c.y), radius_t("t")
{
	gType = TIME_DEPENDENT_CIRCLE ;
	double t1 = c.t - r0/rate ;
	radius_t -= t1 ;
	radius_t *= rate ;
	center.t = t1 ;
}

TimeDependentCircle::TimeDependentCircle(double r0, double rate, Point * c): Circle(r0, c->x, c->y), radius_t("t")
{
	gType = TIME_DEPENDENT_CIRCLE ;
	double t1 = c->t - r0/rate ;
	radius_t -= t1 ;
	radius_t *= rate ;
	center.t = t1 ;
}

double TimeDependentCircle::radiusAtTime(const Point& p) const
{
	return VirtualMachine().eval( radius_t, p) ;
}

Circle TimeDependentCircle::circleAtTime(const Point& p) const
{
	return Circle( radiusAtTime(p), center) ;
}

double TimeDependentCircle::radiusAtTime(Point * p) const
{
	return VirtualMachine().eval( radius_t, p) ;
}

Circle TimeDependentCircle::circleAtTime( Point * p) const
{
	return Circle( radiusAtTime(p), center) ;
}


bool TimeDependentCircle::in(const Point& v) const
{
	return circleAtTime(v).in(v) ;
}

void TimeDependentCircle::project(Point* p) const
{
	return circleAtTime(p).project(p) ;
}

void TimeDependentCircle::setTimeCircles(Vector & instants)
{
	circles.resize(instants.size());
	for(size_t i = 0 ; i < instants.size() ; i++)
	{
		circles[i] = std::make_pair( VirtualMachine().eval( radius_t, Point(0,0,0,instants[i])), instants[i] ) ;
	}
	time_planes = instants.size() ;
}

