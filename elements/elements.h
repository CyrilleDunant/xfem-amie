// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2013
// Author: Ines Jaouadi <ines.jaouadi@epfl.ch>, (C) 2005-2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __ELEMENTS_H_
#define __ELEMENTS_H_

#include <valarray>
#include <vector>
#include <algorithm>
#include <complex>

#include "../geometry/geometry_3D.h"
#include "../geometry/geometry_2D.h"
#include "../physics/physics_base.h"
#include "../polynomial/vm_base.h"
#include "../polynomial/vm_token.h"
#include "../elements/integrable_entity.h"
namespace Amie
{

Function XTransform(const PointArray & points ,const std::valarray< Function> &basis) ;
double xTransform(const Point & p, const PointArray & points, const std::valarray<Function > & basis) ;
Function YTransform(const PointArray & points ,const std::valarray< Function> &basis) ;
double yTransform(const Point & p, const PointArray & points, const std::valarray<Function > & basis) ;
Function ZTransform(const PointArray & points ,const std::valarray< Function> &basis) ;
double zTransform(const Point & p, const PointArray & points, const std::valarray<Function > & basis) ;
Function TTransform(const PointArray & points ,const std::valarray< Function> &basis) ;
double tTransform(const Point & p, const PointArray & points, const std::valarray<Function > & basis) ;
Point coordinateTransform(const Point & p, const PointArray & points, const std::valarray<Function > & basis) ;
Function dXTransform(const PointArray & points ,const std::valarray< Function> &basis, Variable v) ;
Function dYTransform(const PointArray & points ,const std::valarray< Function> &basis, Variable v) ;
Function dZTransform(const PointArray & points ,const std::valarray< Function> &basis, Variable v) ;
Function dTTransform(const PointArray & points ,const std::valarray< Function> &basis, Variable v) ;
double dXTransform(const PointArray & points ,const std::valarray< Function> &basis, Variable v, const Point & p ) ;
double dYTransform(const PointArray & points ,const std::valarray< Function> &basis, Variable v, const Point & p ) ;
double dZTransform(const PointArray & points ,const std::valarray< Function> &basis, Variable v, const Point & p ) ;
double dTTransform(const PointArray & points ,const std::valarray< Function> &basis, Variable v, const Point & p ) ;



class ElementarySurface ;
class NonLinearForm ;

class ElementarySurface : public IntegrableEntity
{

protected:

    std::vector< Function > enrichfunc ;
    std::vector< const Geometry *> enrichmentSource ;
    Form * behaviour ;
    NonLinearForm * nonlinbehaviour ;

public:

    ElementarySurface() ;
    virtual Function jacobian() const = 0;
    virtual ~ElementarySurface() ;
    virtual void clearElementaryMatrix() = 0 ;
    virtual void print()  const = 0 ;

    virtual bool isMoved() const = 0 ;

    virtual const std::valarray< Function > & getShapeFunctions() const = 0;
    virtual std::valarray< Function > & getShapeFunctions() = 0;

    virtual const Geometry * getEnrichmentSource(size_t i) const  {
        return enrichmentSource[i] ;
    } ;

    virtual const std::vector< size_t > getDofIds() const ;

    virtual const Function & getShapeFunction(size_t i) const ;
    virtual Function & getShapeFunction(size_t i) final;

    virtual const Function & getEnrichmentFunction(size_t i) const ;
    virtual const std::vector<Function>  & getEnrichmentFunctions() const final;

    virtual Function getXTransform() const ;
    virtual Function getYTransform() const ;
    virtual Function getTTransform() const final;
    Function getdXTransform(Variable) const ;
    Function getdYTransform(Variable) const ;
    Function getdTTransform(Variable) const ;
    double getdXTransform(Variable, const Point p) const ;
    double getdYTransform(Variable, const Point p) const ;
    double getdTTransform(Variable, const Point p) const ;

    void setEnrichment(const Function & p, const Geometry * g) ;
    virtual Point inLocalCoordinates(const Point & p) const  = 0;

    virtual void getInverseJacobianMatrix(const Point & p, Matrix & ret) = 0 ;

    virtual Vector getNonLinearForces() = 0 ;

    virtual Form * getBehaviour() const final;
    void setBehaviour( Mesh< DelaunayTriangle, DelaunayTreeItem >* msh, Form* f );

    virtual NonLinearForm * getNonLinearBehaviour() const;
    void setNonLinearBehaviour(NonLinearForm * f) ;

    virtual void step(double dt, const Vector * displacements) ;
    virtual void nonLinearStep(double dt, const Vector * displacements) ;

    Order getOrder() const final;
    void setOrder(Order) ;

    virtual void compileAndPrecalculate();
    virtual std::vector<size_t> clearEnrichment(const Geometry * g) final;
    virtual std::vector<size_t> clearAllEnrichment() final;

} ;

class TriElement : public Triangle, public ElementarySurface
{
protected :
    std::valarray<std::valarray<Matrix> > cachedElementaryMatrix ;
    std::valarray<std::valarray<Matrix> > cachedViscousElementaryMatrix ;
    std::valarray< Function > * shapefunc ;
    std::vector<Matrix> cachedJinv ;
    bool isFather ;
protected :

    const GaussPointArray & genGaussPoints();

public:

    bool moved ;

    GEO_DERIVED_OBJECT(Triangle) ;

    TriElement( Point * p0,  Point * p1,  Point * p2) ;

    TriElement(Order order = LINEAR) ;
    void refresh(const TriElement * parent) ;

    virtual std::valarray<std::valarray<Matrix> > & getElementaryMatrix() ;
    virtual std::valarray<std::valarray<Matrix> > & getViscousElementaryMatrix() ;
    virtual void clearElementaryMatrix() final{
        cachedElementaryMatrix.resize(0);
        cachedViscousElementaryMatrix.resize(0) ;
    } ;
    virtual std::valarray<std::valarray<Matrix> > getNonLinearElementaryMatrix(Vector * state)  ;
    virtual void getSecondJacobianMatrix(const Point &p, Matrix & t1, Matrix & t2)  ;
    virtual void getThirdJacobianMatrix(const Point &p, Matrix & t1, Matrix & t2, Matrix & t3) ;

    Function jacobian() const ;

    double  jacobianAtPoint(const Point & p) const ;

    void getInverseJacobianMatrix(const Point & p, Matrix & ret) final;

    const GaussPointArray & getGaussPoints() final;

    virtual const std::valarray< Function > & getShapeFunctions() const final;
    virtual std::valarray< Function > & getShapeFunctions() final;

    virtual bool isMoved() const final;

    virtual void print() const;

    virtual Point inLocalCoordinates(const Point &p) const final;
    virtual std::valarray<std::valarray<Matrix> > getNonLinearElementaryMatrix() ;

    virtual Vector getNonLinearForces()  ;
    virtual Function getXTransform() const final;
    virtual Function getYTransform() const final;
    virtual Function getXTransformAtCentralNodalTime() const final;
    virtual Function getYTransformAtCentralNodalTime() const final;
    virtual Function getZTransformAtCentralNodalTime() const final{
        return Function("0") ;
    }
    virtual Function getTTransformAtCentralNodalTime() const final{
        return Function("0") ;
    }

    virtual ~TriElement();

    virtual Mesh< DelaunayTriangle, DelaunayTreeItem >* get2DMesh() const {
        return nullptr ;
    }
    virtual Mesh< DelaunayTetrahedron, DelaunayTreeItem3D >* get3DMesh() const {
        return nullptr ;
    }

    friend class ElementarySurface ;

} ;

class ElementaryVolume : public IntegrableEntity
{
protected:

    std::vector< Function> enrichfunc ;
    std::vector<Geometry *> enrichmentSource ;
    virtual const GaussPointArray & genGaussPoints()= 0 ;
    Form * behaviour ;
    Order order ;
    NonLinearForm * nonlinbehaviour ;

    std::valarray<std::valarray<Matrix> > cachedElementaryMatrix ;
    std::valarray<std::valarray<Matrix> > cachedViscousElementaryMatrix ;
    Vector cachedForces ;

public:

    ElementaryVolume(bool f = false) ;
    virtual Function jacobian() const ;
    double jacobianAtPoint(const Point & p) const ;
    virtual ~ElementaryVolume() ;

    virtual bool isMoved() const ;

// 	virtual const std::vector<std::pair<size_t,const Function &> > getDofs() const ;
    virtual const std::vector< size_t > getDofIds() const ;
    virtual void clearElementaryMatrix() {
        cachedElementaryMatrix.resize(0);
    } ;
    virtual void print()  const = 0 ;
    virtual const GaussPointArray & getGaussPoints() = 0 ;

    virtual std::valarray<std::valarray<Matrix> > & getElementaryMatrix() = 0;
    virtual Form * getBehaviour() const final;
    virtual void setBehaviour( Mesh< DelaunayTetrahedron, DelaunayTreeItem3D >* msh, Form* f ) final;

    virtual NonLinearForm * getNonLinearBehaviour() const ;

    virtual const std::valarray< Function > & getShapeFunctions() const = 0 ;
    virtual std::valarray< Function > & getShapeFunctions() = 0 ;
    virtual const Function & getShapeFunction(size_t i) const ;
    virtual Function & getShapeFunction(size_t i)  final;
    virtual const Function & getEnrichmentFunction(size_t i) const final ;
// 	virtual Function & getEnrichmentFunction(size_t i)  ;
    virtual const std::vector< Function>  & getEnrichmentFunctions() const final;
// 	virtual std::vector< Function>  & getEnrichmentFunctions() ;
    virtual Function getXTransform() const ;
    virtual Function getYTransform() const ;
    virtual Function getZTransform() const ;
    virtual Function getTTransform() const ;
    virtual Function getdXTransform(Variable v) const ;
    virtual Function getdYTransform(Variable v) const ;
    virtual Function getdZTransform(Variable v) const ;
    virtual Function getdTTransform(Variable v) const ;
    virtual double getdXTransform(Variable v, const Point & p) const ;
    virtual double getdYTransform(Variable v, const Point & p) const ;
    virtual double getdZTransform(Variable v, const Point & p) const ;
    virtual double getdTTransform(Variable v, const Point & p) const ;

    virtual void setEnrichment(const Function &  p, Geometry * g) final;

    virtual void getInverseJacobianMatrix(const Point & p, Matrix & ret) ;

    virtual const PointArray & getBoundingPoints() const = 0;
    virtual PointArray & getBoundingPoints() = 0;
    virtual const Point & getBoundingPoint(size_t i) const = 0;
    virtual Point & getBoundingPoint(size_t i) = 0;
    virtual const PointArray & getInPoints() const = 0;
    virtual PointArray & getInPoints() = 0;
    virtual const Point &  getPoint(size_t i) const = 0 ;
    virtual Point &  getPoint(size_t i) = 0 ;


    virtual void step(double dt, const Vector * displacements) ;
    virtual void nonLinearStep(double dt, const Vector *displacements) ;

    virtual Order getOrder() const final;
    virtual void setOrder(Order) final;

    virtual void compileAndPrecalculate();
    virtual std::vector<size_t> clearEnrichment(const Geometry * g) ;
    virtual std::vector<size_t> clearAllEnrichment() ;

} ;


class TetrahedralElement : public Tetrahedron,  public ElementaryVolume
{
protected :
    std::valarray< Function > *shapefunc ;
    virtual const GaussPointArray & genGaussPoints() final;
    bool isFather ;
    std::vector<Matrix> cachedJinv ;
public:
    bool moved;
    GEO_DERIVED_OBJECT(Tetrahedron) ;

    TetrahedralElement( Point * p0,  Point * p1,  Point * p2, Point * p3) ;
    TetrahedralElement( Point * p0,  Point * p1,  Point * p2, Point * p3, Point * p4,  Point * p5,  Point * p6, Point * p7) ;
    TetrahedralElement(Order order = LINEAR);
    TetrahedralElement(TetrahedralElement * parent, Tetrahedron * t);
    virtual std::valarray<std::valarray<Matrix> > & getElementaryMatrix() ;
    virtual std::valarray<std::valarray<Matrix> > & getViscousElementaryMatrix() ;
    virtual std::valarray<std::valarray<Matrix> > getNonLinearElementaryMatrix() ;
    virtual void getInverseJacobianMatrix(const Point & p, Matrix & ret) ;
    virtual const std::valarray< Function > & getShapeFunctions() const ;
    virtual std::valarray< Function > & getShapeFunctions() ;
    virtual Vector getNonLinearForces() ;

    void refresh(const TetrahedralElement * parent);

    virtual void print() const;
    virtual Point inLocalCoordinates(const Point & p) const final;

    virtual Function getXTransform() const final;
    virtual Function getYTransform() const final;
    virtual Function getZTransform() const final;
    virtual Function getTTransform() const final;


    virtual Mesh< DelaunayTriangle, DelaunayTreeItem >* get2DMesh() const {
        return nullptr ;
    };
    virtual Mesh< DelaunayTetrahedron, DelaunayTreeItem3D >* get3DMesh() const {
        return nullptr ;
    };

    virtual const GaussPointArray & getGaussPoints() final
    {
        return genGaussPoints() ;
    }

    virtual Function getdXTransform(Variable v) const final;
    virtual Function getdYTransform(Variable v) const final;
    virtual Function getdZTransform(Variable v) const final;
    virtual Function getdTTransform(Variable v) const final;

    virtual double getdXTransform(Variable v, const Point & p) const final;
    virtual double getdYTransform(Variable v, const Point & p) const final;
    virtual double getdZTransform(Variable v, const Point & p) const final;
    virtual double getdTTransform(Variable v, const Point & p) const final;

    virtual ~TetrahedralElement();

    friend class ElementaryVolume ;
} ;

class HexahedralElement final: public Hexahedron,  public ElementaryVolume
{
protected :
    const GaussPointArray & genGaussPoints() final;
    std::valarray< Function > *shapefunc ;
    std::valarray<std::valarray<Matrix> > cachedElementaryMatrix ;
public:
    std::vector<HexahedralElement *> neighbourhood ;
// 	bool moved;
    GEO_DERIVED_OBJECT(Hexahedron) ;

    HexahedralElement(Order order, bool f = true) ;
    HexahedralElement(HexahedralElement * parent,Hexahedron * t);

    virtual std::valarray<std::valarray<Matrix> > & getElementaryMatrix() ;
    virtual std::valarray<std::valarray<Matrix> > & getViscousElementaryMatrix() ;


    const GaussPointArray & getGaussPoints() final
    {
        return genGaussPoints() ;
    }
    void refresh(const HexahedralElement * parent);
    virtual void print()  const;
    virtual Point inLocalCoordinates(const Point & p) const ;
    virtual Vector getNonLinearForces() ;
    virtual std::valarray<std::valarray<Matrix> > getNonLinearElementaryMatrix() ;
    virtual const std::valarray< Function > & getShapeFunctions() const ;
    virtual std::valarray< Function > & getShapeFunctions()  ;
    bool visited ;

    virtual Function getXTransform() const ;
    virtual Function getYTransform() const ;
    virtual Function getZTransform() const ;

    virtual Mesh< DelaunayTriangle, DelaunayTreeItem >* get2DMesh() const {
        return nullptr ;
    };
    virtual Mesh< DelaunayTetrahedron, DelaunayTreeItem3D >* get3DMesh() const {
        return nullptr ;
    };

} ;



void computeNeighbourhoodForStructuredHexahedralMesh(std::vector<HexahedralElement *> & vec) ;
void burn(std::vector<HexahedralElement *> & vec) ;
GaussPointArray gaussPointSet(Order order, const TriElement * t) ;
GaussPointArray gaussPointSet(Order order, const TetrahedralElement * t) ;
} 
#endif // __ELEMENTS_H_
