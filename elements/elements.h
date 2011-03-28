// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Ines Jaouadi <ines.jaouadi@epfl.ch>, (C) 2005-2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
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

Mu::Function XTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis) ;
double xTransform(const Mu::Point & p, const Mu::PointArray & points, const std::valarray<Mu::Function > & basis) ;
Mu::Function YTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis) ;
double yTransform(const Mu::Point & p, const Mu::PointArray & points, const std::valarray<Mu::Function > & basis) ;
Mu::Function ZTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis) ;
double zTransform(const Mu::Point & p, const Mu::PointArray & points, const std::valarray<Mu::Function > & basis) ;
Mu::Function TTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis) ;
double tTransform(const Mu::Point & p, const Mu::PointArray & points, const std::valarray<Mu::Function > & basis) ;
Mu::Point coordinateTransform(const Mu::Point & p, const Mu::PointArray & points, const std::valarray<Mu::Function > & basis) ;
Mu::Function dXTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v) ;
Mu::Function dYTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v) ;
Mu::Function dZTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v) ;
Mu::Function dTTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v) ;
double dXTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Mu::Point & p ) ;
double dYTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Mu::Point & p ) ;
double dZTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Mu::Point & p ) ;
double dTTransform(const Mu::PointArray & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Mu::Point & p ) ;

namespace Mu
{

class ElementarySurface ;
class NonLinearForm ;

class ElementarySurface : public IntegrableEntity
{

protected:
	std::valarray< Function > * shapefunc ;
	std::vector< Function > enrichfunc ;
	std::vector< const Geometry *> enrichmentSource ;
	Form * behaviour ;
	NonLinearForm * nonlinbehaviour ;

public:

	ElementarySurface(bool f = false) ;
	virtual Function jacobian() const = 0;
	virtual ~ElementarySurface() ;
	virtual void clearElementaryMatrix() = 0 ;
	virtual void print()  const = 0 ;
	
	bool isFather ;
	
	virtual bool isMoved() const = 0 ;
	
	virtual const std::valarray< Function > & getShapeFunctions() const ;

	virtual const std::vector< size_t > getDofIds() const ;
	
	virtual const Function & getShapeFunction(size_t i) const ;
	
	virtual const Function & getEnrichmentFunction(size_t i) const ;
	virtual const std::vector<Function>  & getEnrichmentFunctions() const ;
	
	virtual Function getXTransform() const ;
	virtual Function getYTransform() const ;
	virtual Function getTTransform() const ;
	Function getdXTransform(Variable) const ;
	Function getdYTransform(Variable) const ;
	Function getdTTransform(Variable) const ;
	double getdXTransform(Variable, const Point p) const ;
	double getdYTransform(Variable, const Point p) const ;
	double getdTTransform(Variable, const Point p) const ;
	
	void setEnrichment(const Function & p, const Geometry * g) ;
	virtual Point inLocalCoordinates(const Point & p) const  = 0;
	
	virtual void getInverseJacobianMatrix(const Point & p, Matrix & ret) const = 0 ;
	
	virtual Vector getNonLinearForces() = 0 ;
		
	virtual Form * getBehaviour() const ;
	void setBehaviour(Form *);
	
	virtual NonLinearForm * getNonLinearBehaviour() const;
	void setNonLinearBehaviour(NonLinearForm * f) ;
	
	virtual void step(double dt, const Vector * displacements) ;
	virtual void stepBack() ;
	virtual void nonLinearStep(double dt, const Vector * displacements) ;
	
	Order getOrder() const;
	void setOrder(Order) ;

	virtual void compileAndPrecalculate();
	virtual std::vector<size_t> clearEnrichment(const Geometry * g) ;
	
} ;


class TriElement : public Triangle, public ElementarySurface
{
protected :
	std::valarray<std::valarray<Matrix> > cachedElementaryMatrix ;
protected :
	
	const GaussPointArray & genGaussPoints();
	
public:
	
	bool moved ;
	
	GEO_DERIVED_OBJECT(Triangle) ;
	
	TriElement( Point * p0,  Point * p1,  Point * p2) ;
	
	TriElement(Order order = LINEAR, bool father = true) ;
	void refresh(const TriElement * parent) ;
	
	virtual std::valarray<std::valarray<Matrix> > & getElementaryMatrix() ;
	virtual void clearElementaryMatrix() { cachedElementaryMatrix.resize(0);} ;
	virtual std::valarray<std::valarray<Matrix> > getNonLinearElementaryMatrix(Vector * state)  ;
	
	Function jacobian() const ;
	
	double  jacobianAtPoint(const Point & p) const ;
	
	void getInverseJacobianMatrix(const Point & p, Matrix & ret) const ;
	
	const GaussPointArray & getGaussPoints();
	
	virtual bool isMoved() const;
	
	virtual void print() const;
	
	virtual Point inLocalCoordinates(const Point &p) const ;
	virtual std::valarray<std::valarray<Matrix> > getNonLinearElementaryMatrix() ;
	
	virtual Vector getNonLinearForces()  ;
	virtual Function getXTransform() const ;
	virtual Function getYTransform() const ;
	
    virtual Mesh< DelaunayTriangle, DelaunayTreeItem >* get2DMesh() const {return NULL ;}
    virtual Mesh< DelaunayTetrahedron, DelaunayTreeItem3D >* get3DMesh() const {return NULL ;}
	
} ;

class ElementaryVolume : public IntegrableEntity
{
protected:
	std::valarray< Function > *shapefunc ;
	std::vector< Function> enrichfunc ;
	std::vector<Geometry *> enrichmentSource ;
	virtual const GaussPointArray & genGaussPoints()= 0 ;
	Form * behaviour ;
	Order order ;
	NonLinearForm * nonlinbehaviour ;
	
	std::valarray<std::valarray<Matrix> > cachedElementaryMatrix ;
	Vector cachedForces ;

public:

	ElementaryVolume(bool f = false) ;
	virtual Function jacobian() const ;
	double jacobianAtPoint(const Point & p) const ;
	virtual ~ElementaryVolume() ;
	
	const bool isFather ;
	
	virtual bool isMoved() const ;
	
// 	virtual const std::vector<std::pair<size_t,const Function &> > getDofs() const ;
	virtual const std::vector< size_t > getDofIds() const ;
	virtual void clearElementaryMatrix() { cachedElementaryMatrix.resize(0);} ;
	virtual void print()  const = 0 ;
	virtual const GaussPointArray & getGaussPoints() = 0 ;

	virtual std::valarray<std::valarray<Matrix> > & getElementaryMatrix() = 0;
	virtual Form * getBehaviour() const ;
	virtual void setBehaviour(Form *);

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
		
	virtual const PointArray & getBoundingPoints() const = 0;
	virtual PointArray & getBoundingPoints() = 0;
	virtual const Point & getBoundingPoint(size_t i) const = 0;
	virtual Point & getBoundingPoint(size_t i) = 0;
	virtual const PointArray & getInPoints() const = 0;
	virtual PointArray & getInPoints() = 0;
	virtual const Point &  getPoint(size_t i) const = 0 ;
	virtual Point &  getPoint(size_t i) = 0 ;
		
	
	virtual void step(double dt, const Vector * displacements) ;
	virtual void stepBack() ;
	virtual void nonLinearStep(double dt, const Vector *displacements) ;
	
	virtual Order getOrder() const;
	virtual void setOrder(Order) ;
	
	virtual void compileAndPrecalculate();
	virtual std::vector<size_t> clearEnrichment(const Geometry * g) ;
} ;


class TetrahedralElement : public Tetrahedron,  public ElementaryVolume
{
protected :
	const GaussPointArray & genGaussPoints();
public:
	bool moved;
	GEO_DERIVED_OBJECT(Tetrahedron) ;
	
	TetrahedralElement( Point * p0,  Point * p1,  Point * p2, Point * p3, bool father = false) ;
	TetrahedralElement( Point * p0,  Point * p1,  Point * p2, Point * p3, Point * p4,  Point * p5,  Point * p6, Point * p7, bool father = false) ;
	TetrahedralElement(Order order = LINEAR, bool father = true);
	TetrahedralElement(TetrahedralElement * parent, Tetrahedron * t);
	virtual std::valarray<std::valarray<Matrix> > & getElementaryMatrix() ;
	virtual std::valarray<std::valarray<Matrix> > getNonLinearElementaryMatrix() ;
	virtual void getInverseJacobianMatrix(const Point & p, Matrix & ret) const;
		
	virtual Vector getNonLinearForces() ;
	
	void refresh(const TetrahedralElement * parent);

	virtual void print() const;
	virtual Point inLocalCoordinates(const Point & p) const ;

	virtual Function getXTransform() const ;
	virtual Function getYTransform() const ;
	virtual Function getZTransform() const ;
	virtual Function getTTransform() const ;

	
	virtual Mesh< DelaunayTriangle, DelaunayTreeItem >* get2DMesh() const {return NULL ;};
	virtual Mesh< DelaunayTetrahedron, DelaunayTreeItem3D >* get3DMesh() const {return NULL ;};
	
	virtual const GaussPointArray & getGaussPoints()
	{
		return genGaussPoints() ;
	}
	
	virtual Function getdXTransform(Variable v) const ;
	virtual Function getdYTransform(Variable v) const ;
	virtual Function getdZTransform(Variable v) const ;
	virtual Function getdTTransform(Variable v) const ;

	virtual double getdXTransform(Variable v, const Point & p) const ;
	virtual double getdYTransform(Variable v, const Point & p) const ;
	virtual double getdZTransform(Variable v, const Point & p) const ;
	virtual double getdTTransform(Variable v, const Point & p) const ;

} ;

class HexahedralElement : public Hexahedron,  public ElementaryVolume
{
protected :
	const GaussPointArray & genGaussPoints() ;
	std::valarray<std::valarray<Matrix> > cachedElementaryMatrix ;
public:
	std::vector<HexahedralElement *> neighbourhood ;
// 	bool moved;
	GEO_DERIVED_OBJECT(Hexahedron) ;
	
	HexahedralElement(Order order, bool f = true) ;
	HexahedralElement(HexahedralElement * parent,Hexahedron * t);

	virtual std::valarray<std::valarray<Matrix> > & getElementaryMatrix() ;
	
	const GaussPointArray & getGaussPoints()
	{
		return genGaussPoints() ;
	}
	void refresh(const HexahedralElement * parent);
	virtual void print()  const;
	virtual Point inLocalCoordinates(const Point & p) const ;
	virtual Vector getNonLinearForces() ;
	virtual std::valarray<std::valarray<Mu::Matrix> > getNonLinearElementaryMatrix() ;
	bool visited ;

	virtual Function getXTransform() const ;
	virtual Function getYTransform() const ;
	virtual Function getZTransform() const ;
	
	virtual Mesh< DelaunayTriangle, DelaunayTreeItem >* get2DMesh() const {return NULL ;};
    virtual Mesh< DelaunayTetrahedron, DelaunayTreeItem3D >* get3DMesh() const {return NULL ;};

} ;

} ;

void computeNeighbourhoodForStructuredHexahedralMesh(std::vector<Mu::HexahedralElement *> & vec) ;
void burn(std::vector<Mu::HexahedralElement *> & vec) ;

#endif // __ELEMENTS_H_
