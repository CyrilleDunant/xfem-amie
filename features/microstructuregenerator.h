
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
#include "../utilities/random.h"

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

/** \brief Utility class to convert circles to any distribution*/
struct InclusionConverter
{
    GeometryType geom ;
    RandomDistribution * area ;
    RandomDistribution * aspectRatio ;
    RandomDistribution * orientation ;

    InclusionConverter(GeometryType type, RandomDistribution * a = new ConstantDistribution(1.), RandomDistribution * ar = new ConstantDistribution(1.), RandomDistribution * o = new ConstantDistribution(0.)) ;

    void setArea(RandomDistribution * a) ;
    void setAspectRatio(RandomDistribution * ar) ;
    void setOrientation(RandomDistribution * o) ;

    void setArea(double a) ;
    void setAspectRatio(double ar) ;
    void setOrientation(double o) ;

    Feature * convert(Inclusion * inc) const ;
    std::vector<Feature *> convert(std::vector<Inclusion *> inc) const ;


} ;

struct InclusionGenerator
{
    double authorizeRotationsDuringPlacement ;

    InclusionGenerator(double rot = 0.) : authorizeRotationsDuringPlacement(rot) { } ;

    virtual Feature * convert(Inclusion * inc) const { return inc ; }
    virtual std::vector<Feature *> convert(std::vector<Inclusion *> inc) const ;

} ;

struct EllipsoidalInclusionGenerator : public InclusionGenerator
{
    double shape ;
    double orientation ;
    double orientationVariability ;
    double shapeVariability ;

    EllipsoidalInclusionGenerator( double s, double o = 0., double ov = M_PI , double sv = 0., double rot = 0.) : InclusionGenerator(rot), shape(s), orientation(o), orientationVariability(ov), shapeVariability(sv) { } ;

    virtual Feature * convert(Inclusion * inc) const ;

} ;

struct RectangularInclusionGenerator : public InclusionGenerator
{
    double shape ;
    double orientation ;
    double orientationVariability ;
    double shapeVariability ;

    RectangularInclusionGenerator( double s, double o = 0., double ov = M_PI , double sv = 0., double rot = 0.) : InclusionGenerator(rot), orientation(o), orientationVariability(ov), shapeVariability(sv) { } ;

    virtual Feature * convert(Inclusion * inc) const ;

} ;

struct PolygonalInclusionGenerator : public InclusionGenerator
{
    int vertex ;
    double orientation ;
    double orientationVariability ;
    int vertexVariability ;

    PolygonalInclusionGenerator(int v, double o = 0., double ov = M_PI, int vv = 0, double rot = 0.) : InclusionGenerator(rot), vertex(v), orientation(o), orientationVariability(ov), vertexVariability(vv) { } ;

    virtual PolygonalSample * generatePolygon(double radius, size_t npoints, double phase = 0) const;
    virtual Feature * convert(Inclusion * inc ) const ;
} ;

// see Beddow and Meloy 1980, cited by Wang et al 1999
struct GravelPolygonalInclusionGenerator : public PolygonalInclusionGenerator
{
    double p ;
    double b ;
    size_t m ;

    GravelPolygonalInclusionGenerator(double p_, double b_, size_t m_, int v, double o = 0., double ov = M_PI, int vv = 0, double rot = 0.) : PolygonalInclusionGenerator(v,o,ov,vv,rot), p(p_), b(b_),m(m_) { } ;

    virtual PolygonalSample * generatePolygon(double radius, size_t npoints, double phase = 0) const;
} ;

// see Wang et al 1999 option A
struct CrushedPolygonalInclusionGenerator : public PolygonalInclusionGenerator
{
    double shape ; // between 0 and 1

    CrushedPolygonalInclusionGenerator(double s, int v, double o = 0., double ov = M_PI, int vv = 0, double rot = 0.) : PolygonalInclusionGenerator(v,o,ov,vv, rot), shape(s<0? -s : s) { if(s > 1-POINT_TOLERANCE) { shape = 0.999 ; } } ;

    virtual PolygonalSample * generatePolygon(double radius, size_t npoints, double phase = 0) const;
} ;

// see Wang et al 1999 option B
struct CrushedSubtendedPolygonalInclusionGenerator : public PolygonalInclusionGenerator
{
    double shape ; // between 0 and 1
    double delta ; // between 0 and 1

    CrushedSubtendedPolygonalInclusionGenerator(double s, double d, int v, double o = 0., double ov = M_PI, int vv = 0, double rot = 0.) : PolygonalInclusionGenerator(v,o,ov,vv,rot), shape(s<0? -s : s), delta(d<0? -d: d) { if(shape > 1-POINT_TOLERANCE) { shape = 0.999 ; } if(delta > 1-POINT_TOLERANCE) { delta = 0.999 ; } } ;

    virtual PolygonalSample * generatePolygon(double radius, size_t npoints, double phase = 0) const ;
} ;

}





#endif // MICROSTRUCTURE_GENERATOR_H
