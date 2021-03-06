//
// C++ Interface: inclusion3d
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2006-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __INCLUSION3D_H__
#define __INCLUSION3D_H__

#include "features.h"

namespace Amie
{

/** \brief 3D spherical inclusion*/
class Inclusion3D :  public Sphere,  public Feature
{
    friend class Sphere ;
public:

    /** \brief Construct inclusion from radius and center position
    *
    * @param father Father feature
    * @param radius of the inclusion
    * @param x center x
    * @param y center y
    * @param z center z
    */
    Inclusion3D(Feature *father, double radius, double x, double y, double z) ;

    /** \brief Construct inclusion from radius and center position
    *
    * @param father Father feature
    * @param radius of the inclusion
    * @param center center
    */
    Inclusion3D(Feature *father, double radius,  Point center) ;

    /** \brief Construct inclusion from radius and center position
    *
    * @param radius of the inclusion
    * @param x center x
    * @param y center y
    * @param z center z
    */
    Inclusion3D(double radius, double x, double y, double z) ;

    /** \brief Construct inclusion from radius and center position
    *
    * @param radius of the inclusion
    * @param center center
    */
    Inclusion3D(double radius, Point center) ;

    /** \brief do nothing */
    virtual void addSamplePoints(PointSet * po ) { };

    /** \brief return true if the boundary overlaps that of the argument*/
    virtual bool interacts(Feature * f, double d) const ;

    /** \brief get list of refinement zones*/
    virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;

    /** \brief return empty vector*/
    virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  ;

    /** \brief return all tets in mesh with at least a vertex in this Feature*/
    virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt)  ;

    virtual void print() const
    {
        std::cout << "I am an inclusion" << std::endl ;
    }

    virtual void setRadius(double newr);

    /** \brief return false*/
    virtual bool isVoid( const Point &) const {
        return false ;
    }


public:

    GEO_DERIVED_OBJECT(Sphere) ;

    virtual void sample(double linearDensity, double surfaceDensityFactor, Sampler * sampler = nullptr) ;

} ;


/** \brief Regular octahedron inclusion*/
class OctahedralInclusion :  public RegularOctahedron,  public Feature
{
    friend class RegularOctahedron ;
public:

    /** \brief Construct inclusion from length and center position
    *
    * @param father Father feature
    * @param length of the inclusion
    * @param x center x
    * @param y center y
    * @param z center z
    */
    OctahedralInclusion(Feature *father, double radius, double x, double y, double z) ;

    /** \brief Construct inclusion from length and center position
    *
    * @param father Father feature
    * @param length of the inclusion
    * @param center center
    */
    OctahedralInclusion(Feature *father, double radius,  Point center) ;

    /** \brief Construct inclusion from length and center position
    *
    * @param length of the inclusion
    * @param x center x
    * @param y center y
    * @param z center z
    */
    OctahedralInclusion(double radius, double x, double y, double z) ;

    /** \brief Construct inclusion from length size and center position
    *
    * @param length of the inclusion
    * @param center center
    */
    OctahedralInclusion(double length, Point center) ;

    /** \brief do nothing */
    virtual void addSamplePoints(PointSet * po ) { };

    /** \brief return true if the boundary overlaps that of the argument*/
    virtual bool interacts(Feature * f, double d) const ;

    /** \brief get list of refinement zones*/
    virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;

    /** \brief return empty vector*/
    virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  ;

    /** \brief return all tets in mesh with at least a vertex in this Feature*/
    virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt)  ;

    virtual void print() const
    {
        std::cout << "I am an inclusion" << std::endl ;
    }

    /** \brief return false*/
    virtual bool isVoid( const Point &) const {
        return false ;
    }


public:

    GEO_DERIVED_OBJECT(RegularOctahedron) ;

    virtual void sample(double linearDensity, double surfaceDensityFactor, Sampler * sampler = nullptr) ;

} ;

/** \brief Inclusion of spherical shape used only to attribute behaviour. This feature has no effect on the mesh*/
class VirtualInclusion3D :  public Sphere,  public VirtualFeature
{
    friend class Sphere ;
public:

    /** \brief Construct inclusion from radius and center position
    *
    * @param father Father feature
    * @param radius of the inclusion
    * @param x center x
    * @param y center y
    * @param z center z
    */
    VirtualInclusion3D(Feature *father, double radius, double x, double y, double z) ;

    /** \brief Construct inclusion from radius and center position
    *
    * @param father Father feature
    * @param radius of the inclusion
    * @param center center
    */
    VirtualInclusion3D(Feature *father, double radius,  Point center) ;

    /** \brief Construct inclusion from radius and center position
    *
    * @param radius of the inclusion
    * @param x center x
    * @param y center y
    * @param z center z
    */
    VirtualInclusion3D(double radius, double x, double y, double z) ;

    /** \brief Construct inclusion from radius and center position
    *
    * @param radius of the inclusion
    * @param center center
    */
    VirtualInclusion3D(double radius, Point center) ;

    /** \brief do nothing*/
    virtual void addSamplePoints(PointSet * po ) { };

    /** \brief return false*/
    virtual bool interacts(Feature * f, double d) const ;

    /** \brief return empty vector*/
    virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;

    /** \brief return empty vector*/
    virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  ;

    /** \brief return all tets in mesh with at least a vertex in this Feature*/
    virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt)  ;

    /** \brief return tthis*/
    virtual Feature * getSource() ;

    /** \brief return false */
    virtual bool inBoundary(const Point *, double d) const ;

    /** \brief return false*/
    virtual bool inBoundary(const Point &, double d) const ;

    virtual void print() const
    {
        std::cout << "I am an inclusion" << std::endl ;
    }

    /** \brief return false*/
    virtual bool isVoid( const Point &) const {
        return false ;
    }


public:

    GEO_DERIVED_OBJECT(Sphere) ;

    virtual void sample(double linearDensity, double surfaceDensityFactor, Sampler * sampler = nullptr) ;

} ;

}

#endif
