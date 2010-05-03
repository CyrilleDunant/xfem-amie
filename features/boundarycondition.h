// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010

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

public:
	BoundaryCondition(LagrangeMultiplierType t, const double & d) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const = 0 ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t) const = 0 ;
	void setData(double newval) { data = newval ;}
	double getData() const { return data ;}
} ;

/** \brief Boundary condition object for usage in multigrid solver. Work in Progress*/
class ProjectionDefinedBoundaryCondition : public BoundaryCondition
{
private:
	Point direction ;

public:
	ProjectionDefinedBoundaryCondition(LagrangeMultiplierType t,const Point & direction, double d = 0) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  const ;
} ;

/** \brief Boundary condition object for usage in multigrid solver.*/
class BoundingBoxDefinedBoundaryCondition : public BoundaryCondition
{
private:
	BoundingBoxPosition pos ;
	
public:
	BoundingBoxDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double d = 0 ) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  const ;
} ;

/** \brief Boundary condition object for usage in multigrid solver.*/
class BoundingBoxAndRestrictionDefinedBoundaryCondition : public BoundaryCondition
{
private:
	BoundingBoxPosition pos ;
	double xmin, xmax, ymin, ymax, zmin, zmax ;
	
public:
	BoundingBoxAndRestrictionDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp,double  ym, double yp, double zm, double zp, double d = 0 ) ;
	BoundingBoxAndRestrictionDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp,double  ym, double yp, double d = 0 ) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  const ;
} ;

class BoundingBoxNearestNodeDefinedBoundaryCondition : public BoundaryCondition
{
private:
	BoundingBoxPosition pos ;
	Point nearest ;
	
public:
	BoundingBoxNearestNodeDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, Point nearest, double d = 0 ) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  const ;
} ;

/** \brief Boundary condition object for usage in multigrid solver*/
class GeometryDefinedBoundaryCondition : public BoundaryCondition
{
private:
	Geometry * domain ;

public:
	GeometryDefinedBoundaryCondition(LagrangeMultiplierType t, Geometry * source, double d = 0) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  const ;
} ;

/** \brief Boundary condition object for usage in multigrid solver. Work in Progress*/
class GeometryProjectedBoundaryCondition : public BoundaryCondition
{
private:
	Geometry * domain ;
	Point from ;
	Point direction ;
public:
	GeometryProjectedBoundaryCondition(LagrangeMultiplierType t, Geometry * source, const Point & from,  const Point & direction, double d = 0 ) ;
	virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const ;
	virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  const ;
} ;
} ;

#endif
