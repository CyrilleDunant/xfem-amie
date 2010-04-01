// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
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
#include "../elements/integrable_entity.h"

Mu::Function XTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis) ;
double xTransform(const Mu::Point & p, const std::valarray<Mu::Point*> & points, const std::valarray<Mu::Function > & basis) ;
Mu::Function YTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis) ;
double yTransform(const Mu::Point & p, const std::valarray<Mu::Point*> & points, const std::valarray<Mu::Function > & basis) ;
Mu::Function ZTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis) ;
double zTransform(const Mu::Point & p, const std::valarray<Mu::Point*> & points, const std::valarray<Mu::Function > & basis) ;
Mu::Function TTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis) ;
double tTransform(const Mu::Point & p, const std::valarray<Mu::Point*> & points, const std::valarray<Mu::Function > & basis) ;
Mu::Point coordinateTransform(const Mu::Point & p, const std::valarray<Mu::Point*> & points, const std::valarray<Mu::Function > & basis) ;
Mu::Function dXTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v) ;
Mu::Function dYTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v) ;
Mu::Function dZTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v) ;
Mu::Function dTTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v) ;
double dXTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Mu::Point & p ) ;
double dYTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Mu::Point & p ) ;
double dZTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Mu::Point & p ) ;
double dTTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Mu::Point & p ) ;

namespace Mu
{

class ElementarySurface ;
class NonLinearForm ;

class ElementarySurface : public IntegrableEntity
{
protected:
	std::valarray< Function > * shapefunc ;
	std::vector< Function > enrichfunc ;
	std::vector< Geometry *> enrichmentSource ;
	Form * behaviour ;
	NonLinearForm * nonlinbehaviour ;
	std::vector<std::vector<Matrix> > cachedElementaryMatrix ;
	GaussPointArray cachedGaussPoints ;
	Vector cachedForces ;

public:

	ElementarySurface(bool f = false) ;
	virtual Function jacobian() const = 0;
	virtual ~ElementarySurface() ;
	
	virtual void print()  const = 0 ;
	
	bool isFather ;
	
	virtual bool isMoved() const = 0 ;
	
	virtual const std::valarray< Function > & getShapeFunctions() const ;

// 	virtual const std::vector<std::pair<size_t,const Function &> > getDofs() const ;
	virtual const std::vector< size_t > getDofIds() const ;
	
	virtual const Function & getShapeFunction(size_t i) const ;
// 	virtual Function & getShapeFunction(size_t i)  ;
	
	virtual const Function & getEnrichmentFunction(size_t i) const ;
// 	virtual Function & getEnrichmentFunction(size_t i) ;
	virtual const std::vector<Function>  & getEnrichmentFunctions() const ;
// 	virtual std::vector<Function>  & getEnrichmentFunctions() ;
	
	virtual Function getXTransform() const ;
	virtual Function getYTransform() const ;
	virtual Function getTTransform() const ;
	Function getdXTransform(Variable) const ;
	Function getdYTransform(Variable) const ;
	Function getdTTransform(Variable) const ;
	double getdXTransform(Variable, const Point p) const ;
	double getdYTransform(Variable, const Point p) const ;
	double getdTTransform(Variable, const Point p) const ;
	
	void setEnrichment(const Function & p, Geometry * g) ;
	virtual Point inLocalCoordinates(const Point & p) const  = 0;
	
	virtual void getInverseJacobianMatrix(const Point & p, Matrix & ret) const = 0 ;
	
	virtual std::vector<std::vector<Matrix> > & getElementaryMatrix() = 0;
	virtual Vector getForces() { return Vector(0) ;}
	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix()  = 0;
	virtual Vector getNonLinearForces() = 0 ;
	
	virtual const Point &  getPoint(size_t i) const = 0 ;
	virtual  Point &  getPoint(size_t i)  = 0 ;
	
	virtual Form * getBehaviour() const ;
	void setBehaviour(Form *);
	
	virtual NonLinearForm * getNonLinearBehaviour() const;
	void setNonLinearBehaviour(NonLinearForm * f) ;
	
	void step(double dt, Vector * displacements) ;
	void stepBack() ;
	void nonLinearStep(double dt, Vector *displacements) ;
	
	Order getOrder() const;
	void setOrder(Order) ;

	virtual void compileAndPrecalculate();
	virtual void clearEnrichment(const Geometry * g) ;
	
} ;


class TriElement : public Triangle, public ElementarySurface
{
	
protected :
	
	GaussPointArray genGaussPoints() const;
	virtual void computeCenter();
	
public:
	
	bool moved ;
	
	GEO_DERIVED_OBJECT(Triangle) ;
	
	TriElement( Point * p0,  Point * p1,  Point * p2, bool father = false) ;
	
	TriElement(Order order = LINEAR, bool father = true) ;
	void refresh(const TriElement * parent) ;
	
	virtual std::vector<std::vector<Matrix> > & getElementaryMatrix() ;
	
	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix(Vector * state)  ;
	
	Function jacobian() const ;
	
	double  jacobianAtPoint(const Point p) const ;
	
	void getInverseJacobianMatrix(const Point & p, Matrix & ret) const ;
	
	GaussPointArray getGaussPoints() const;
	
	virtual bool isMoved() const;
	
	virtual void print() const;
	
	virtual Point inLocalCoordinates(const Point &p) const ;
	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix() ;
	
	virtual Vector getNonLinearForces()  ;
	virtual Function getXTransform() const ;
	virtual Function getYTransform() const ;
	
} ;

// class QuadElement : public ConvexGeometry, public ElementarySurface
// {
// protected :
// 	
// 	GaussPointArray genGaussPoints() const;
// 	virtual void computeCenter();
// 	
// public:
// 	
// 	bool moved ;
// 	
// 	GEO_DERIVED_OBJECT(ConvexGeometry) ;
// 	
// 	QuadElement( Point * p0,  Point * p1,  Point * p2,  Point * p3, bool father = false) ;
// 	
// 	QuadElement(Order order = LINEAR, bool father = true) ;
// 	void refresh(const QuadElement * parent) ;
// 	
// 	virtual std::vector<std::vector<Matrix> > & getElementaryMatrix() ;
// 	
// 	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix(Vector * state)  ;
// 	
// 	Function jacobian() const ;
// 	
// 	double  jacobianAtPoint(const Point p) const ;
// 	
// 	void getInverseJacobianMatrix(const Point & p, Matrix & ret) const ;
// 	
// 	GaussPointArray getGaussPoints() const;
// 	
// 	virtual bool isMoved() const;
// 	
// 	virtual void print() const;
// 	
// 	virtual Point inLocalCoordinates(const Point &p) const ;
// 	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix() ;
// 	
// 	virtual Vector getNonLinearForces()  ;
// 	virtual Function getXTransform() const ;
// 	virtual Function getYTransform() const ;
// } ;

class ElementaryVolume : public IntegrableEntity
{
protected:
	std::valarray< Function > *shapefunc ;
	std::vector< Function> enrichfunc ;
	std::vector<Geometry *> enrichmentSource ;
	virtual GaussPointArray genGaussPoints() const = 0 ;
	Form * behaviour ;
	Order order ;
	NonLinearForm * nonlinbehaviour ;
	
	std::vector<std::vector<Matrix> > cachedElementaryMatrix ;
	Vector cachedForces ;
	GaussPointArray cachedGaussPoints ;

public:

	ElementaryVolume(bool f = false) ;
	virtual Function jacobian() const ;
	double jacobianAtPoint(const Point & p) const ;
	virtual ~ElementaryVolume() ;
	
	const bool isFather ;
	
	virtual bool isMoved() const ;
	
// 	virtual const std::vector<std::pair<size_t,const Function &> > getDofs() const ;
	virtual const std::vector< size_t > getDofIds() const ;
	virtual void clearElementaryMatrix() { } ;
	virtual void print()  const = 0 ;
	virtual GaussPointArray getGaussPoints() const { return genGaussPoints() ;}
	virtual std::vector<std::vector<Matrix> > & getElementaryMatrix() = 0;
	virtual Form * getBehaviour() const ;
	virtual void setBehaviour(Form *);
	virtual Vector getForces() { return Vector(0) ;}

	virtual NonLinearForm * getNonLinearBehaviour() const ;
	
	virtual const std::valarray< Function > & getShapeFunctions() const ;
	virtual const Function & getShapeFunction(size_t i) const ;
// 	virtual Function & getShapeFunction(size_t i)  ;
	virtual const Function & getEnrichmentFunction(size_t i) const  ;
// 	virtual Function & getEnrichmentFunction(size_t i)  ;
	virtual const std::vector< Function>  & getEnrichmentFunctions() const ;
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

	virtual void setEnrichment(const Function &  p, Geometry * g) ;
	
	virtual void getInverseJacobianMatrix(const Point & p, Matrix & ret) const ;
		
	virtual const std::valarray< Point * > & getBoundingPoints() const = 0;
	virtual std::valarray< Point * > & getBoundingPoints() = 0;
	virtual const Point & getBoundingPoint(size_t i) const = 0;
	virtual Point & getBoundingPoint(size_t i) = 0;
	virtual const std::valarray< Point * > & getInPoints() const = 0;
	virtual std::valarray< Point * > & getInPoints() = 0;
	virtual const Point &  getPoint(size_t i) const = 0 ;
	virtual Point &  getPoint(size_t i) = 0 ;
		
	
	virtual void step(double dt, Vector * displacements) ;
	virtual void stepBack() ;
	virtual void nonLinearStep(double dt, Vector *displacements) ;
	
	virtual Order getOrder() const;
	virtual void setOrder(Order) ;
	
	virtual void compileAndPrecalculate();
	virtual void clearEnrichment(const Geometry * g) ;
} ;


class TetrahedralElement : public Tetrahedron,  public ElementaryVolume
{
protected :
	GaussPointArray genGaussPoints() const;
	virtual void computeCenter();
public:
	bool moved;
	GEO_DERIVED_OBJECT(Tetrahedron) ;
	
	TetrahedralElement( Point * p0,  Point * p1,  Point * p2, Point * p3, bool father = false) ;
	TetrahedralElement( Point * p0,  Point * p1,  Point * p2, Point * p3, Point * p4,  Point * p5,  Point * p6, Point * p7, bool father = false) ;
	TetrahedralElement(Order order = LINEAR, bool father = true);
	TetrahedralElement(TetrahedralElement * parent, Tetrahedron * t);
	virtual std::vector<std::vector<Matrix> > & getElementaryMatrix() ;
	
	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix() ;
	
	virtual Vector getForces() ;
	
	virtual Vector getNonLinearForces() ;
	
	void refresh(const TetrahedralElement * parent);

	virtual void print() const;
	virtual Point inLocalCoordinates(const Point & p) const ;

	virtual Function getXTransform() const ;
	virtual Function getYTransform() const ;
	virtual Function getZTransform() const ;

} ;

class HexahedralElement : public Hexahedron,  public ElementaryVolume
{
protected :
	GaussPointArray genGaussPoints() const;
	std::vector<std::vector<Matrix> > cachedElementaryMatrix ;
	virtual void computeCenter();
public:
	std::vector<HexahedralElement *> neighbourhood ;
// 	bool moved;
	GEO_DERIVED_OBJECT(Hexahedron) ;
	
	HexahedralElement(Order order, bool f = true) ;
	HexahedralElement(HexahedralElement * parent,Hexahedron * t);

	virtual std::vector<std::vector<Matrix> > & getElementaryMatrix() ;

	virtual Vector getForces() const ;

	void refresh(const HexahedralElement * parent);
	virtual void print()  const;
	virtual Point inLocalCoordinates(const Point & p) const ;
	virtual Vector getNonLinearForces() ;
	virtual std::vector<std::vector<Mu::Matrix> > getNonLinearElementaryMatrix() ;
	bool visited ;

	virtual Function getXTransform() const ;
	virtual Function getYTransform() const ;
	virtual Function getZTransform() const ;

} ;

} ;

void computeNeighbourhoodForStructuredHexahedralMesh(std::vector<Mu::HexahedralElement *> & vec) ;
void burn(std::vector<Mu::HexahedralElement *> & vec) ;

#endif // __ELEMENTS_H_
