// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011

#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include "../geometry/geometry_base.h"
#include "../elements/elements.h" 
#include "../mesher/delaunay.h"
#include "../mesher/delaunay_3d.h"
#include "../mesher/structuredmesh.h"
#include "../utilities/samplingcriterion.h"
#include "../utilities/grid.h"
#include "../solvers/assembly.h"
#include "feature_base.h"

namespace Mu {
/** \brief Abstract boundary condition object for usage in multigrid solver.
 * 
*/
class BoundaryCondition
{	
protected:
	LagrangeMultiplierType condition;
	double data ;
	double scale ;
	Function dataFunction ;
	int axis ;
	bool function ;
	std::vector<DelaunayTriangle *> cache2d ;
	std::vector<DelaunayTetrahedron *> cache3d ;
	std::vector<std::vector<Point> > cache ;
public:
	BoundaryCondition(LagrangeMultiplierType t, const double & d, int a = 0) ;
	BoundaryCondition(LagrangeMultiplierType t, const Function & d, int a = 0) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) = 0 ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t) = 0 ;
	void setData(double newval) { data = newval ;}
	void setData(const Function & f) { dataFunction = f ;}
	double getData() const { return data ;}
	int getAxisIndex() const { return axis ; }
	void setAxisIndex(int a) { axis = a ; }
	void clearCache()
	{
		cache.clear();
		cache2d.clear();
		cache3d.clear();
	} ;
	
	virtual void setScale(double) ;
	virtual double getScale() const;
		
	LagrangeMultiplierType getConditionType() const {return condition ;}  ;
	virtual ~BoundaryCondition() {} ;
} ;

class PlaneSectionsBoundaryConditions : public BoundaryCondition
{
	bool isVertical ;
	double uplimit ;
	double downlimit ;
public:
	PlaneSectionsBoundaryConditions(bool isVertical, double down, double up, int a = 0) :BoundaryCondition(GENERAL, 0.),  isVertical(isVertical), downlimit(down), uplimit(up) { };
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t) ;
} ;


class NullBoundaryCondition : public BoundaryCondition
{
	public:
		NullBoundaryCondition() : BoundaryCondition(nullptr_CONDITION, 0.) { } ;
		virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) {} ;
		virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t) {} ;
} ;

/** \brief Boundary condition object for usage in multigrid solver. Work in Progress*/
class ProjectionDefinedBoundaryCondition : public BoundaryCondition
{
private:
	Point direction ;

public:
	ProjectionDefinedBoundaryCondition(LagrangeMultiplierType t,const Point & direction, double d = 0, int a = 0) ;
	ProjectionDefinedBoundaryCondition(LagrangeMultiplierType t,const Point & direction, const Function & d, int a = 0 ) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

/** \brief Boundary condition object for usage in multigrid solver.*/
class BoundingBoxDefinedBoundaryCondition : public BoundaryCondition
{
private:
	BoundingBoxPosition pos ;
	
public:
	BoundingBoxDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double d = 0, int a = 0 ) ;
	BoundingBoxDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, const Function & d, int a = 0 ) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

/** \brief Boundary condition object for usage in multigrid solver.*/
class BoundingBoxAndRestrictionDefinedBoundaryCondition : public BoundaryCondition
{
private:
	BoundingBoxPosition pos ;
	double xmin, xmax, ymin, ymax, zmin, zmax ;
	
public:
	BoundingBoxAndRestrictionDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp,double  ym, double yp, double zm, double zp, double d = 0, int a = 0 ) ;
	BoundingBoxAndRestrictionDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp,double  ym, double yp, double d = 0, int a = 0 ) ;
	BoundingBoxAndRestrictionDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp,double  ym, double yp, double zm, double zp, const Function & d, int a = 0 ) ;
	BoundingBoxAndRestrictionDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp,double  ym, double yp, const Function & d, int a = 0 ) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

class BoundingBoxNearestNodeDefinedBoundaryCondition : public BoundaryCondition
{
private:
	BoundingBoxPosition pos ;
	Point nearest ;
	
public:
	BoundingBoxNearestNodeDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, Point nearest, double d = 0, int a = 0 ) ;
	BoundingBoxNearestNodeDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, Point nearest, const Function & d, int a = 0 ) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

/** \brief Boundary condition applied on points on the surface of specified geometry*/
class GeometryDefinedBoundaryCondition : public BoundaryCondition
{
protected:
	Geometry * domain ;

public:
	GeometryDefinedBoundaryCondition(LagrangeMultiplierType t, Geometry * source, double d = 0, int a = 0) ;
	GeometryDefinedBoundaryCondition(LagrangeMultiplierType t, Geometry * source, const Function & d, int a = 0) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

/** \brief Boundary condition applied on points in specified geometry*/
class GeometryDefinedSurfaceBoundaryCondition : public BoundaryCondition
{
protected:
	Geometry * domain ;
public:
	GeometryDefinedSurfaceBoundaryCondition(LagrangeMultiplierType t, Geometry * source, double d = 0, int a = 0) ;
	GeometryDefinedSurfaceBoundaryCondition(LagrangeMultiplierType t, Geometry * source, const Function & d, int a = 0) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

/** \brief Boundary condition applied on points in specified geometry*/
class GeometryDefinedSurfaceAndRestrictionsBoundaryCondition : public BoundaryCondition
{
protected:
	Geometry * domain ;
	double xmin, xmax, ymin, ymax, zmin, zmax ;
public:
	GeometryDefinedSurfaceAndRestrictionsBoundaryCondition(LagrangeMultiplierType t, Geometry * source, double d = 0, int a = 0) ;
	GeometryDefinedSurfaceAndRestrictionsBoundaryCondition(LagrangeMultiplierType t, Geometry * source, const Function & d, int a = 0) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

/** \brief Boundary condition object for usage in multigrid solver. Work in Progress*/
class GeometryProjectedBoundaryCondition : public BoundaryCondition
{
private:
	Geometry * domain ;
	Point from ;
	Point direction ;
public:
	GeometryProjectedBoundaryCondition(LagrangeMultiplierType t, Geometry * source, const Point & from,  const Point & direction, double d = 0 , int a = 0) ;
	GeometryProjectedBoundaryCondition(LagrangeMultiplierType t, Geometry * source, const Point & from,  const Point & direction, const Function & d , int a = 0) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;


/** \brief Boundary condition for MAD: take the displacement field on the border of an element, and apply it to the mesh.*/
class ElementDefinedBoundaryCondition : public BoundaryCondition
{
	protected:
		ElementarySurface * surface ;
		ElementaryVolume * volume ;
	public:
		ElementDefinedBoundaryCondition(ElementarySurface * surface) ;
		ElementDefinedBoundaryCondition(ElementaryVolume * volume) ;
		virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
		virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

/** \brief Boundary condition for single Dofs*/
class DofDefinedBoundaryCondition : public BoundaryCondition
{
	protected:
		size_t id ;
		ElementarySurface * surface ;
		ElementaryVolume * volume ;
		std::valarray<Matrix> * Jinv ;
		GaussPointArray * gp ;
	public:
		DofDefinedBoundaryCondition(LagrangeMultiplierType t, ElementarySurface * surface, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, size_t id, double d = 0, int a = 0 ) ;
		DofDefinedBoundaryCondition(LagrangeMultiplierType t, ElementaryVolume * surface, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv , size_t id, double d = 0, int a = 0 ) ;
		DofDefinedBoundaryCondition(LagrangeMultiplierType t, ElementarySurface * surface, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv , size_t id, const Function & d, int a = 0 ) ;
		DofDefinedBoundaryCondition(LagrangeMultiplierType t, ElementaryVolume * surface, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv , size_t id, const Function & d, int a = 0 ) ;
		virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
		virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
		virtual ~DofDefinedBoundaryCondition() { delete gp ; delete Jinv ;}
} ;

class TimeContinuityBoundaryCondition : public BoundaryCondition
{
	protected:
		std::vector<std::pair<size_t, size_t> > correspondance ;
  
	public:
		Vector previousDisp ;
		bool goToNext ;
		double initialValue ;
		TimeContinuityBoundaryCondition(double initialValue = 0) ;
		virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
		virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
		
} ;

class GlobalForceBoundaryCondition : public BoundaryCondition
{
	protected:
		Vector dataVector ;

	public:
		GlobalForceBoundaryCondition(Vector & data) ;

		void setDataVector(Vector & d) ;
		virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
		virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;


} ;

#endif
