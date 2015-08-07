// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#ifndef __ST_GEOMETRY_2D_H_
#define __ST_GEOMETRY_2D_H_

#include "geometry_2D.h"
#include "../polynomial/vm_function_base.h"

namespace Amie
{

class TimeDependentCircle : public Circle
{
protected:
    Function radius_t ;
    std::vector< std::pair<double, double> > circles ;
    bool linear ;
    bool constant ;

public:
    TimeDependentCircle(Function r = Function(), const Point & center = Point());
    TimeDependentCircle(Function r, Point* center);
    TimeDependentCircle(double r0, double rate, const Point & center);
    TimeDependentCircle(double r0, double rate, Point * center);

    void setTimeCircles( Vector & instants ) ;
    void setLinear(bool l) {
        linear = l ;
        constant = false ;
    }
    void setConstant(bool c) {
        linear = false ;
        constant = c ;
    }

    double radiusAtTime(const Point & p) const ;
    double radiusAtTime(Point * p) const ;
    Circle circleAtTime(const Point & p) const ;
    Circle circleAtTime(Point * p) const ;

    double timeAtRadius( double r ) const ;

    virtual bool in(const Point &v) const ;

    virtual void setRadius(double ) { };

    virtual void project(Point * p) const;

    Function getRadiusFunction() const {
        return radius_t ;
    }
    Function getRadiusFunction(Function & time) const ;

    void setInitialTime( double t) {
        center.getT() = t ;
    }

    void setRadiusFunction(Function & f) { radius_t = f ; }

} ;

}

#endif // _ST_GEOMETRY_2D_H_
