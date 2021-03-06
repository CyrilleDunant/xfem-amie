#include "space_time_geometry_2D.h"
#include "../polynomial/vm_base.h"
#include <iomanip>

namespace Amie{

TimeDependentCircle::TimeDependentCircle(Function r, const Point & c): Circle(0., c.getX(), c.getY()), radius_t(r)
{
	gType = TIME_DEPENDENT_CIRCLE ;
	time_planes = 2 ;
	center.getT() = c.getT() ;
	linear = false ;
	constant = false ;
}


TimeDependentCircle::TimeDependentCircle(Function r, Point * c): Circle(0., c->getX(), c->getY()), radius_t(r)
{
	gType = TIME_DEPENDENT_CIRCLE ;
	time_planes = 2 ;
	center.getT() = c->getT() ;
	linear = false ;
	constant = false ;
}

TimeDependentCircle::TimeDependentCircle(double r0, double rate, const Point & c): Circle(r0, c.getX(), c.getY()), radius_t("t")
{
	gType = TIME_DEPENDENT_CIRCLE ;
	time_planes = 2 ;
	double t1 = c.getT() - r0/rate ;
	radius_t -= t1 ;
	radius_t *= rate ;
	center.getT() = t1 ;
	linear = true ;
	constant = false ;
}

TimeDependentCircle::TimeDependentCircle(double r0, double rate, Point * c): Circle(r0, c->getX(), c->getY()), radius_t("t")
{
	gType = TIME_DEPENDENT_CIRCLE ;
	time_planes = 2 ;
	double t1 = c->getT() - r0/rate ;
	radius_t -= t1 ;
	radius_t *= rate ;
	center.getT() = t1 ;
	linear = true ;
	constant = false ;
}

double TimeDependentCircle::radiusAtTime(const Point& p) const
{
// 	if(VirtualMachine().eval( radius_t, p) == 0)
// 		p.print() ;
	return VirtualMachine().eval( radius_t, p) ;
}

Circle TimeDependentCircle::circleAtTime(const Point& p) const
{
	Point c = center ;
	c.getT() = p.getT() ;
	return Circle( radiusAtTime(p), c) ;
}

double TimeDependentCircle::radiusAtTime(Point * p) const
{
	return VirtualMachine().eval( radius_t, p) ;
}

Circle TimeDependentCircle::circleAtTime( Point * p) const
{
	Point c = center ;
	c.getT() = p->getT() ;
	return Circle( radiusAtTime(p), center) ;
}


bool TimeDependentCircle::in(const Point& v) const
{
	if(v.getT() < center.getT())
		return false ;
	return circleAtTime(v).in(v) ;
}

void TimeDependentCircle::project(Point* p) const
{
	if(p->getT() < center.getT()) 
	{
		p->getX() = center.getX() ;
		p->getY() = center.getY() ;
		return ;
	}

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

Function TimeDependentCircle::getRadiusFunction(Function& time) const
{
	Function r = radius_t ;
// 	r.setNumberOfVariableTransforms(4);
 	r.setVariableTransform(TIME_VARIABLE, time);
	r.makeVariableTransformDerivative() ;
//	r.setNumberOfDerivatives(0);
	return r ;
}

double TimeDependentCircle::timeAtRadius(double r) const
{
	if(constant)
		return -1 ;
	if(linear)
	{
		double r0 =  VirtualMachine().eval( radius_t, 0,0,0,0) ;
		double r1 =  VirtualMachine().eval( radius_t, 0,0,0,1) ;
		return (r-r0)/(r1-r0) ;
	}
	else
	{
		double r0 =  VirtualMachine().eval( radius_t, 0,0,0,0) ;
		double r1 =  VirtualMachine().eval( radius_t, 0,0,0,1) ;
		return (r-r0)/(r1-r0) ;
	}
}

}

