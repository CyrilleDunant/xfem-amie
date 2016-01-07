
/*
    mesh abstract implementation for AMIE
    Copyright (C) 2010  Cyrille Dunant

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

#ifndef MICROSTRUCTURE_GENERATOR_H
#define MICROSTRUCTURE_GENERATOR_H

#include "feature_base.h"
#include "inclusion.h"
#include "polygonSample.h"
#include "../geometry/geometry_base.h"

namespace Amie
{
struct MicrostructureGenerator
{

    MicrostructureGenerator() { };
    virtual std::vector<Feature *> getFeatures(Geometry * sample) = 0;
    virtual double score() = 0;
    virtual void print() const = 0 ;
} ;


struct AggregateDistribution2DGenerator : public MicrostructureGenerator
{

    double area ;
    double dmax ;
    double fill ;
    double minMaxRatio ;
    double itzSize ;
    double inclusionNumber ;
    double massOfAggregates ;
    Geometry * sample ;

    AggregateDistribution2DGenerator(double area, double dmax, double itzSize, double fill, double minMaxRatio) ;

    virtual std::vector<Feature *> getFeatures(Geometry * sample) ;
    virtual double score() ;
    virtual void print() const ;
} ;

struct InclusionGenerator
{
    std::default_random_engine rng ;
    double authorizeRotationsDuringPlacement ;

    InclusionGenerator(double rot = 0.) : authorizeRotationsDuringPlacement(rot) { } ;

    virtual Feature * convert(Inclusion * inc)  { return inc ; }
    virtual std::vector<Feature *> convert(std::vector<Inclusion *> inc)  ;

} ;

/*PARSE Circular InclusionGenerator
    @value[placement_rotation] 0 // maximum angle of the rotation authorized during placement (if 0, the inclusion will keep its generated orientation)
*/
struct CircularInclusionGenerator : public InclusionGenerator
{
    CircularInclusionGenerator(double rot = 0.) : InclusionGenerator(rot) { } ;
} ;

/*PARSE XFEM InclusionGenerator */
struct XFEMInclusionGenerator : public InclusionGenerator
{
    XFEMInclusionGenerator() : InclusionGenerator(0.) { } ;

    virtual Feature * convert(Inclusion * inc)  ;

} ;

/*PARSE SpaceTimeXFEM InclusionGenerator */
struct SpaceTimeXFEMInclusionGenerator : public XFEMInclusionGenerator
{
    SpaceTimeXFEMInclusionGenerator() : XFEMInclusionGenerator() { } ;

    virtual Feature * convert(Inclusion * inc)  ;

} ;

/*PARSE Ellipsoidal InclusionGenerator
    @value[shape_factor] // ratio between major and minor axis of the ellipses
    @value[orientation] 0 // default angle of the major axis
    @value[orientation_variability] 3.141592 // extent of the random uniform distribution of the angles of the major axis
    @value[shape_factor_variability] 0 // extent of the random uniform distribution of the shape factor
    @value[placement_rotation] 0 // maximum angle of the rotation authorized during placement (if 0, the inclusion will keep its generated orientation)
*/
struct EllipsoidalInclusionGenerator : public InclusionGenerator
{
    std::uniform_real_distribution< double > shape ;
    std::uniform_real_distribution< double > orientation ;

    EllipsoidalInclusionGenerator( double s, double o = 0., double ov = M_PI , double sv = 0., double rot = 0.) : InclusionGenerator(rot), shape(s-sv,s+sv), orientation(o-ov, o+ov) { } ;

    virtual Feature * convert(Inclusion * inc)  ;

} ;

/*PARSE Rectangular InclusionGenerator
    @value[shape_factor] // ratio between major and minor sides of the rectangle
    @value[orientation] 0 // default angle of the major axis
    @value[orientation_variability] 3.141592 // extent of the random uniform distribution of the angles of the major axis
    @value[shape_factor_variability] 0 // extent of the random uniform distribution of the shape factor
    @value[placement_rotation] 0 // maximum angle of the rotation authorized during placement (if 0, the inclusion will keep its generated orientation)
*/
struct RectangularInclusionGenerator : public InclusionGenerator
{
    std::uniform_real_distribution< double > shape ;
    std::uniform_real_distribution< double > orientation ;

    RectangularInclusionGenerator( double s, double o = 0., double ov = M_PI , double sv = 0., double rot = 0.) : InclusionGenerator(rot), shape(s-sv,s+sv), orientation(o-ov, o+ov) { } ;

    virtual Feature * convert(Inclusion * inc) ;

} ;

/*PARSE Polygonal InclusionGenerator
    @value[vertex] // number of vertexes of the polygons
    @value[orientation] 0 // default angle of the major axis
    @value[orientation_variability] 3.141592 // extent of the random uniform distribution of the angles of the major axis
    @value[vertex_variability] 0 // extent of the random uniform distribution of the number of vertexes
    @value[placement_rotation] 0 // maximum angle of the rotation authorized during placement (if 0, the inclusion will keep its generated orientation)
*/
struct PolygonalInclusionGenerator : public InclusionGenerator
{
    std::uniform_int_distribution< int > vertex ;
    std::uniform_real_distribution< double > orientation ;
    bool forceOrientation ;

    PolygonalInclusionGenerator(int v, double o = 0., double ov = M_PI, int vv = 0, double rot = 0., bool force = true) : InclusionGenerator(rot), vertex(std::max(v-vv,3),std::max(v+vv,3)), orientation(o-ov,o+ov), forceOrientation(force) { } ;

    virtual PolygonalSample * generatePolygon(double radius) ;
    virtual Feature * convert(Inclusion * inc )  ;
} ;

/*PARSE VoronoiPolygonal InclusionGenerator
    @value[box_width] // width of the box in which the Voronoi distribution is generated
    @value[grains] // number of points from which the Voronoi distribution is generated
    @value[spacing] // distance between two different points of the Voronoi distribution
    @value[shape_factor] 1 // elongation of the particles
    @value[orientation] 0 // default angle of the major axis
    @value[orientation_variability] 3.141592 // extent of the random uniform distribution of the angles of the major axis
    @value[shape_factor_variability] 0 // extent of the random uniform distribution of the shape factor
    @value[placement_rotation] 0 // maximum angle of the rotation authorized during placement (if 0, the inclusion will keep its generated orientation)
*/
struct VoronoiPolygonalInclusionGenerator : public PolygonalInclusionGenerator
{
    std::vector<PolygonalSample *> source ;
    std::uniform_int_distribution< size_t > index ;

    VoronoiPolygonalInclusionGenerator( double box, size_t seed, double minDist, double shape = 1., double o = 0., double ov = M_PI, double sv = 0., double rot = 0., bool force = true) ;

    virtual PolygonalSample * generatePolygon(double radius) ;
} ; 

// see Beddow and Meloy 1980, cited by Wang et al 1999
/*PARSE GravelPolygonal InclusionGenerator
    @value[amplitude_factor] 0.9 // controls the shape of the aggregates
    @value[amplitude_exponent] 1.9 // controls the shape of the aggregates
    @value[degree] 2 // degree of the harmonic function
    @value[vertex] // number of vertexes of the polygons
    @value[orientation] 0 // default angle of the major axis
    @value[orientation_variability] 3.141592 // extent of the random uniform distribution of the angles of the major axis
    @value[vertex_variability] 0 // extent of the random uniform distribution of the number of vertexes
    @value[placement_rotation] 0 // maximum angle of the rotation authorized during placement (if 0, the inclusion will keep its generated orientation)
*/
struct GravelPolygonalInclusionGenerator : public PolygonalInclusionGenerator
{
    double p ;
    double b ;
    size_t m ;
    std::uniform_real_distribution< double > phi ;

    GravelPolygonalInclusionGenerator(double p_, double b_, size_t m_, int v, double o = 0., double ov = M_PI, int vv = 0, double rot = 0., bool f = true) : PolygonalInclusionGenerator(v,o,ov,vv,rot, f), p(p_), b(b_),m(m_), phi(0, 2.*M_PI) { } ;

    virtual PolygonalSample * generatePolygon(double radius) ;
} ;

// see Wang et al 1999 option A
/*PARSE CrushedPolygonal InclusionGenerator
    @value[shape_factor] // ratio between major and minor axis of the aggregate
    @value[vertex] // number of vertexes of the polygons
    @value[orientation] 0 // default angle of the major axis
    @value[orientation_variability] 3.141592 // extent of the random uniform distribution of the angles of the major axis
    @value[vertex_variability] 0 // extent of the random uniform distribution of the number of vertexes
    @value[placement_rotation] 0 // maximum angle of the rotation authorized during placement (if 0, the inclusion will keep its generated orientation)
*/
struct CrushedPolygonalInclusionGenerator : public PolygonalInclusionGenerator
{
    std::uniform_real_distribution< double > shape ;
    std::uniform_real_distribution< double > phi ;

    CrushedPolygonalInclusionGenerator(double s, int v, double o = 0., double ov = M_PI, int vv = 0, double rot = 0., bool f = true) : PolygonalInclusionGenerator(v,o,ov,vv, rot, f), shape(-(1.-s)/(1.+s),(1.-s)/(1.+s)), phi(0,2.*M_PI) { } ;

    virtual PolygonalSample * generatePolygon(double radius) ;
} ;

// see Wang et al 1999 option B
/*PARSE CrushedSubtendedPolygonal InclusionGenerator
    @value[shape_factor] // ratio between major and minor axis of the aggregate
    @value[angle_variability] // extent of the random uniform distribution of the angles between consecutive segments of the aggregate
    @value[vertex] // number of vertexes of the polygons
    @value[orientation] 0 // default angle of the major axis
    @value[orientation_variability] 3.141592 // extent of the random uniform distribution of the angles of the major axis
    @value[vertex_variability] 0 // extent of the random uniform distribution of the number of vertexes
    @value[placement_rotation] 0 // maximum angle of the rotation authorized during placement (if 0, the inclusion will keep its generated orientation)
*/
struct CrushedSubtendedPolygonalInclusionGenerator : public PolygonalInclusionGenerator
{
    std::uniform_real_distribution< double > shape ;
    std::uniform_real_distribution< double > delta ;

    CrushedSubtendedPolygonalInclusionGenerator(double s, double d, int v, double o = 0., double ov = M_PI, int vv = 0, double rot = 0., bool f = true) : PolygonalInclusionGenerator(v,o,ov,vv,rot, f), shape(-(1.-s)/(1.+s),(1.-s)/(1.+s)), delta(-d, d) { } ;

    virtual PolygonalSample * generatePolygon(double radius)  ;
} ;

}





#endif // MICROSTRUCTURE_GENERATOR_H
