//
// C++ Interface: sample
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2006-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __POLY_SAMPLE_3D_H__
#define __POLY_SAMPLE_3D_H__

#include "features.h"

namespace Amie
{

class PolygonalSample3D :  public PolygonPrism,  public Feature
{
public:

    PolygonalSample3D(Amie::Feature* father, const std::valarray< Amie::Point* >& points, const Amie::Point& vector,const Amie::Point& origin) ;

    virtual bool interacts(Feature * f, double d) const ;

    /** \brief return empty list*/
    virtual std::vector<Geometry *> getRefinementZones(size_t) const {
        return std::vector<Geometry *>() ;
    }

    /** \brief return all triangles in mesh with at least a vertex in this Feature*/
    virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  {
        return std::vector<DelaunayTriangle *>(0) ;
    }

    /** \brief return empty vector*/
    virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt) ;

    virtual void print() const
    {
        std::cout << "I am a polygonal sample" << std::endl ;
    }

    /** \brief return false */
    virtual bool isVoid( const Point &p) const {
        return false;
    }

public:

    GEO_DERIVED_OBJECT(PolygonPrism) ;

    virtual void sample(double linearDensity, double surfaceDensityFactor, Sampler * sampler = nullptr)
    {
        sampleSurface(linearDensity, surfaceDensityFactor, sampler) ;
    }

} ;


} 

#endif
