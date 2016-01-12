// Copyright: See COPYING file that comes with this distribution

#ifndef __SAMPLER_H_
#define __SAMPLER_H_
#include "../geometry_base.h"
#include "../geometry_2D.h"
#include "../geometry_3D.h"

namespace Amie
{

/*SOURCE Sampler */

/*PARSE . Sampler */
class Sampler
{
public:
    Sampler() { } ;

    virtual std::vector<Point> sampleBoundingSurface( const Geometry * geom, double linearDensity ) ;
    virtual std::vector<Point> sampleInnerSurface( const Geometry * geom, double linearDensity ) ;

    virtual std::vector<Point> sampleCircleBoundingSurface( const Circle * geom, double linearDensity ) { return std::vector<Point>() ; }
    virtual std::vector<Point> sampleCircleInnerSurface( const Circle * geom, double linearDensity ) { return std::vector<Point>() ; }

    virtual std::vector<Point> sampleEllipseBoundingSurface( const Ellipse * geom, double linearDensity ) { return std::vector<Point>() ; } 
    virtual std::vector<Point> sampleEllipseInnerSurface( const Ellipse * geom, double linearDensity ) { return std::vector<Point>() ; } 

    virtual std::vector<Point> sampleTriangleBoundingSurface( const Triangle * geom, double linearDensity ) { return std::vector<Point>() ; } 
    virtual std::vector<Point> sampleTriangleInnerSurface( const Triangle * geom, double linearDensity ) { return std::vector<Point>() ; } 

    virtual std::vector<Point> sampleRectangleBoundingSurface( const Rectangle * geom, double linearDensity ) { return std::vector<Point>() ; } 
    virtual std::vector<Point> sampleRectangleInnerSurface( const Rectangle * geom, double linearDensity ) { return std::vector<Point>() ; } 

    virtual std::vector<Point> samplePolygonBoundingSurface( const Polygon * geom, double linearDensity ) { return std::vector<Point>() ; } 
    virtual std::vector<Point> samplePolygonInnerSurface( const Polygon * geom, double linearDensity ) { return std::vector<Point>() ; } 

    virtual std::vector<Point> sampleSphereBoundingSurface( const Sphere * geom, double linearDensity ) { return std::vector<Point>() ; } 
    virtual std::vector<Point> sampleSphereInnerSurface( const Sphere * geom, double linearDensity ) { return std::vector<Point>() ; } 

    virtual std::vector<Point> sampleHexahedronBoundingSurface( const Hexahedron * geom, double linearDensity ) { return std::vector<Point>() ; } 
    virtual std::vector<Point> sampleHexahedronInnerSurface( const Hexahedron * geom, double linearDensity ) { return std::vector<Point>() ; } 

    virtual std::vector<Point> samplePolygonPrismBoundingSurface( const PolygonPrism * geom, double linearDensity ) { return std::vector<Point>() ; } 
    virtual std::vector<Point> samplePolygonPrismInnerSurface( const PolygonPrism * geom, double linearDensity ) { return std::vector<Point>() ; } 

    virtual std::vector<Point> sampleLoftedPolygonBoundingSurface( const LoftedPolygonPrism * geom, double linearDensity ) { return std::vector<Point>() ; } 
    virtual std::vector<Point> sampleLoftedPolygonInnerSurface( const LoftedPolygonPrism * geom, double linearDensity ) { return std::vector<Point>() ; } 

} ;

}

#endif // __SAMPLER_H_
