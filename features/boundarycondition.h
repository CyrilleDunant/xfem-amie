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

namespace Amie {
/** \brief Abstract boundary condition object for usage in multigrid solver.
 *
*/

struct LinearInterpolatedMaterialLaw ;
class GeometryBasedEffect ;
class CollisionDetector ;

class BoundaryCondition
{
    friend Form ;
protected:
    LagrangeMultiplierType condition;
    double data ;
    double scale ;
    bool active ;
    Function dataFunction ;
    LinearInterpolatedMaterialLaw * dataInterpolation = nullptr ;
    int axis ;
    bool function ;
    std::vector<DelaunayTriangle *> cache2d ;
    std::vector<DelaunayTetrahedron *> cache3d ;
    std::vector<std::vector<Point> > cache ;
    
    CollisionDetector * collisionDetection ;
    GeometryBasedEffect * contactLaw  ;

public:
    BoundaryCondition(LagrangeMultiplierType t, const double & d, int a = 0) ;
    BoundaryCondition(LagrangeMultiplierType t, const Function & d, int a = 0) ;
    BoundaryCondition(LagrangeMultiplierType t, CollisionDetector * contactCondition, GeometryBasedEffect * contactLaw ) ;
    
    virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t);
    
    virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t);
    
    virtual void setData(double newval) {
        data = newval ;
    }
    virtual void setData(const Function & f) {
        dataFunction = f ;
    }
    virtual double getData() const {
        return data ;
    }
    
    virtual void setInterpolation( LinearInterpolatedMaterialLaw * inter ) ;

    virtual void step(double relatime, double dt) ;

    virtual void setInterpolation( std::string f ) ; 

    virtual const Function & getDataFunction() const {
        return dataFunction ;
    }
    virtual int getAxisIndex() const {
        return axis ;
    }
    virtual void setAxisIndex(int a) {
        axis = a ;
    }
    virtual void clearCache()
    {
        cache.clear();
        cache2d.clear();
        cache3d.clear();
    }
    virtual void setActive(bool a) {active = a ;}

    virtual void setScale(double) ;
    virtual double getScale() const;

    LagrangeMultiplierType getConditionType() const {
        return condition ;
    } 
    virtual ~BoundaryCondition() {}
} ;

class PlaneSectionsBoundaryConditions final: public BoundaryCondition
{
    bool isVertical ;
    double uplimit ;
    double downlimit ;
public:
    PlaneSectionsBoundaryConditions(bool isVertical, double down, double up, int a = 0) :BoundaryCondition(GENERAL, 0.),  isVertical(isVertical), uplimit(up), downlimit(down) { }
    virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t) ;
} ;


class NullBoundaryCondition final: public BoundaryCondition
{
public:
    NullBoundaryCondition() : BoundaryCondition(nullptr_CONDITION, 0.) { }
    virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) {}
    virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t) {}
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
class BoundingBoxDefinedBoundaryCondition final: public BoundaryCondition
{
private:
    BoundingBoxPosition pos ;

public:
    BoundingBoxDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double d = 0, int a = 0 ) ;
    BoundingBoxDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, CollisionDetector * contactCondition, GeometryBasedEffect * contactLaw ) ;
    BoundingBoxDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, const Function & d, int a = 0 ) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

typedef enum {
    LOADING,
    UNLOADING,
    ULTIMATE_STRAIN,
    ULTIMATE_STRESS
} LoadingState;

class LoadingCycle
{
protected :
    LoadingState type ;
    FeatureTree * ft ;
    double ultimate ;
    double condition ;
    size_t axisIndex ;
    double currentState ;
    double rate ;
    bool cycleAtEnd ;
    bool cycleStarted ;
    LoadingCycle * chainedCycle ;
public:
    LoadingCycle(FeatureTree *ft, LoadingState condition, double ultimate, size_t axisIndex, double initialState, double rate = 0.01) :type(LOADING), ft(ft), ultimate(ultimate), condition(condition), axisIndex(axisIndex), currentState(initialState), rate(rate), cycleAtEnd(false), cycleStarted(false), chainedCycle(nullptr){}
    
    LoadingCycle(FeatureTree *ft, LoadingState condition, double ultimate, size_t axisIndex, LoadingCycle * c, double rate = 0.01) :type(LOADING), ft(ft), ultimate(ultimate), condition(condition), axisIndex(axisIndex), currentState(0), rate(rate), cycleAtEnd(false), cycleStarted(false), chainedCycle(c) {}
    
    double getValue() ;
    double isAtEnd() const ;
} ;

class BoundingBoxCycleDefinedBoundaryCondition final: public BoundaryCondition
{
private:
    std::vector<BoundingBoxPosition> positions ;
    std::vector<LagrangeMultiplierType> types ;
    std::vector<LoadingCycle> cycles ;
    int currentCycle ;
    BoundaryCondition * currentBC ;

public:
    BoundingBoxCycleDefinedBoundaryCondition(std::vector<LoadingCycle> cycles, const std::vector<LagrangeMultiplierType> t, const std::vector<BoundingBoxPosition> & pos) ;
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
    BoundingBoxAndRestrictionDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double zm, double zp,CollisionDetector * contactCondition, GeometryBasedEffect * contactLaw )  ;
    virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

class BoundingBoxNearestNodeDefinedBoundaryCondition final: public BoundaryCondition
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
class GeometryDefinedSurfaceBoundaryCondition final: public BoundaryCondition
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
class GeometryAndFaceDefinedSurfaceBoundaryCondition final: public BoundaryCondition
{
protected:
    Geometry * domain ;
    Point faceNormal ;
public:
    GeometryAndFaceDefinedSurfaceBoundaryCondition(LagrangeMultiplierType t, Geometry * source, const Point & normal, double d = 0, int a = 0) ;
    GeometryAndFaceDefinedSurfaceBoundaryCondition(LagrangeMultiplierType t, Geometry * source, const Point & normal, const Function & d, int a = 0) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

/** \brief Boundary condition applied on points in specified geometry*/
class GeometryDefinedSurfaceAndRestrictionsBoundaryCondition final: public BoundaryCondition
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
class GeometryProjectedBoundaryCondition final: public BoundaryCondition
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
class ElementDefinedBoundaryCondition final: public BoundaryCondition
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
class DofDefinedBoundaryCondition final: public BoundaryCondition
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
    virtual ~DofDefinedBoundaryCondition() {
        delete gp ;
        delete Jinv ;
    }
} ;

class TimeContinuityBoundaryCondition final: public BoundaryCondition
{
protected:
    std::vector<std::pair<size_t, size_t> > correspondance ;

public:
    Vector previousDisp ;
    bool goToNext ;
    double initialValue ;
    double instant ;
    double minDeltaTime ;
    TimeContinuityBoundaryCondition(double initialValue = 0) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;

} ;

class GlobalForceBoundaryCondition final: public BoundaryCondition
{
protected:
    Vector dataVector ;

public:
    GlobalForceBoundaryCondition(Vector & data) ;

    void setDataVector(Vector & d) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;

class GlobalBoundaryCondition final: public BoundaryCondition
{

public:
    
    GlobalBoundaryCondition(LagrangeMultiplierType t, double d = 0, int a = 0) ;
    GlobalBoundaryCondition(LagrangeMultiplierType t,  const Function & d, int a = 0) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) ;
    virtual void apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  ;
} ;


}

#endif
