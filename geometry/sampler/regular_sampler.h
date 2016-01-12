// Copyright: See COPYING file that comes with this distribution

#ifndef __REGULAR_SAMPLER_H_
#define __REGULAR_SAMPLER_H_
#include "sampler.h"

namespace Amie
{

class RegularSampler : public Sampler
{
protected:
    bool force ;

public:
    RegularSampler(bool f = false) : Sampler(), force(f) { } ;

    virtual std::vector<Point> sampleInnerSurface( const Geometry * geom, double linearDensity ) ;

    virtual std::vector<Point> sampleRectangleBoundingSurface( const Rectangle * geom, double linearDensity ) ;
    virtual std::vector<Point> sampleRectangleInnerSurface( const Rectangle * geom, double linearDensity ) ;

} ;

}

#endif // __REGULAR_SAMPLER_H_
