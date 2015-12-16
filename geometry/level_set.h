// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#ifndef __LEVEL_SET_H_
#define __LEVEL_SET_H_
#include "../utilities/matrixops.h"
#include "../polynomial/vm_function_base.h"
#include "geometry_base.h"

#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif
#ifdef HAVE_SSE4
#include <smmintrin.h>
#endif

#include <map>

namespace Amie
{

class LevelSet : public Geometry
{
  
  protected:
    
	std::vector<Function> distanceFunction ;
	Geometry * source ;
	PointArray boundingPoints ;
	
  public:
  
	LevelSet() ;
	LevelSet(Point p) ;
	LevelSet(double x, double y) ;
	LevelSet(Geometry * g) ;
	LevelSet(Function f) ;
	LevelSet(std::vector<Function> f) ;

	virtual ~LevelSet() { } ;

	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_TWO_DIMENSIONAL ;
	}
/*	virtual std::vector<Point> getBoundingBox() const { return std::vector<Point>(0) ; } ;*/
	virtual const PointArray & getBoundingPoints() const { return boundingPoints ; } ;
	virtual PointArray & getBoundingPoints() { return boundingPoints ; } ;
	virtual const Point & getBoundingPoint(size_t i) const {return center ; } ;
	virtual Point & getBoundingPoint(size_t i)  {return center ; } ;
	virtual const PointArray & getInPoints() const {return boundingPoints ; } ;
	virtual PointArray & getInPoints() {return boundingPoints ; } ;
	virtual const Point & getInPoint(size_t i) const {return center ; } ;
	virtual Point & getInPoint(size_t i) {return center ; } ;
	virtual void setBoundingPoint(size_t i, Point * p) { } ;
	virtual void setBoundingPoints(const PointArray & nb) { } ;
	virtual void setInPoints(const PointArray & nb) { } ;
	virtual const Point & getPoint(size_t i) const {return center ; } ;
	virtual Point & getPoint(size_t i)  {return center ; } ;
	virtual Geometry * getSource() {return source ; } ;
	virtual const Point & getCenter() const {return center ; };
	virtual Point & getCenter() {return center ; } ;
	virtual void project(Point *) const { } ;
	virtual void setCenter(const Point & newCenter) {this->center = newCenter ; } ;
	virtual GeometryType getGeometryType() const {return LEVEL_SET ; } ;
	virtual double getRadius() const {return 0 ; } ;
	virtual void sampleBoundingSurface(size_t num_points) { } ;
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const {return std::vector<Point>(0) ; } ;
	virtual void sampleSurface(size_t num_points) { } ;
	virtual bool in(const Point & p)const {return false ; } ;
	virtual size_t size() const {return 0 ; } ;
	virtual double area() const {return 0 ; } ;
	virtual double volume() const {return 0 ; } ;
	virtual size_t sides() const {return 0 ; } ;
	virtual bool intersects(const Geometry *) const {return false ; } ;
	virtual std::vector<Point> intersection(const Geometry *) const {return std::vector<Point>(0) ; } ;
	virtual std::vector<Point> getBoundingBox() const ;
//	virtual size_t timePlanes() const {return 1 ; } ;
//	virtual size_t & timePlanes() {return 1 ; } ;

} ;

}

#endif // __LEVEL_SET_H_
