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
#ifndef __POLY_SAMPLE_H__
#define __POLY_SAMPLE_H__

#include "features.h"

namespace Amie
{

class PolygonalSample :  public Polygon,  public Feature
{
public:

    PolygonalSample(Feature *father, const std::valarray<Point *> & points) ;

    virtual bool interacts(Feature * f, double d) const ;

    /** \brief return empty list*/
    virtual std::vector<Geometry *> getRefinementZones(size_t) const {
        return std::vector<Geometry *>() ;
    }

    /** \brief return all triangles in mesh with at least a vertex in this Feature*/
    virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  ;

    /** \brief return empty vector*/
    virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt) {
        return std::vector<DelaunayTetrahedron *>(0) ;
    }

    virtual void print() const
    {
        std::cout << "I am a polygonal sample" << std::endl ;
    }

    /** \brief return false */
    virtual bool isVoid( const Point &p) const {
        return false;
    }

public:

    GEO_DERIVED_OBJECT(Polygon) ;

    virtual void sample(size_t n)
    {
        this->sampleSurface(n) ;
    }

} ;


}

#endif
