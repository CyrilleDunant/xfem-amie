// Copyright: See COPYING file that comes with this distribution

#ifndef __GRADIENT_SAMPLER_H_
#define __GRADIENT_SAMPLER_H_
#include "sampler.h"

namespace Amie
{

class GradientSampler : public Sampler
{
    Segment direction ;
    double start ;
    double end ;
    bool valid = false ;

public:
    GradientSampler(Segment dir, double s, double e = 1.) : Sampler(), direction(dir), start(s), end(e) { valid = (direction.norm() > POINT_TOLERANCE) ; } ;
    GradientSampler(Point A, Point B, double s, double e = 1.) : Sampler(), direction(A,B), start(s), end(e) { valid = (direction.norm() > POINT_TOLERANCE) ; } ;

    virtual std::vector<Point> sampleRectangleBoundingSurface( const Rectangle * geom, double linearDensity ) ;
    virtual std::vector<Point> sampleRectangleInnerSurface( const Rectangle * geom, double linearDensity ) ;

    double getLinearDensity( Point p ) ;

} ;

}

#endif // __GRADIENT_SAMPLER_H_
