// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution

#ifndef __GEOMETRY_WITH_EXCLUSION_H_
#define __GEOMETRY_WITH_EXCLUSION_H_
#include "geometry_base.h"
#include "../utilities/matrixops.h"

namespace Amie
{
class GeometryWithExclusion : public Geometry
{
protected:
	Geometry * father ;
	
public:
	std::vector<Geometry *> exclusions ;

	GeometryWithExclusion(Geometry * f, Geometry * e) ;
	GeometryWithExclusion(Geometry * f, std::vector<Geometry *> e) ;
	
	virtual void project(Point *) const ;
	virtual bool in(const Point & p) const  ;
	virtual double area() const ;
	virtual double volume() const ;
	virtual bool intersects(const Geometry *) const ;
	virtual std::vector<Point> intersection(const Geometry *) const ;
	virtual void computeCenter() { }
	virtual const PointArray& getBoundingPoints() const { return father->getBoundingPoints() ; }
        virtual PointArray& getBoundingPoints() { return father->getBoundingPoints() ; }
        virtual const Point& getBoundingPoint(size_t i) const { return father->getBoundingPoint(i) ; }
        virtual Point& getBoundingPoint(size_t i)  { return father->getBoundingPoint(i) ; }
        virtual void setBoundingPoint(size_t i, Point* p)  { return father->setBoundingPoint(i, p) ; }
        virtual void setBoundingPoints(const PointArray& p) { return father->setBoundingPoints(p) ; }
        virtual const Point& getPoint(size_t i) const  { return father->getPoint(i) ; }
        virtual Point& getPoint(size_t i)  { return father->getPoint(i) ; }
        virtual double getRadius() const { return father->getRadius() ; }
        virtual void sampleBoundingSurface(size_t n) { }
        virtual std::vector<Point> getSamplingBoundingPoints(size_t) const { std::vector<Point> ret ; return ret ; }
        virtual void sampleSurface(size_t) { }
        virtual size_t size() const { return father->size() ; } 
        virtual std::vector<Point, std::allocator<Point> > getBoundingBox() const { return father->getBoundingBox() ; }
        virtual SpaceDimensionality spaceDimensions() const { return father->spaceDimensions() ; }
	
} ;

} ;

#endif
