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
#include "../physics/physics.h"
#include "../polynomial/vm_base.h"
#include "../elements/integrable_entity.h"

Mu::Function XTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis) ;
Mu::Function YTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis) ;
Mu::Function ZTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis) ;
Mu::Function TTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis) ;
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
	std::vector< std::pair< size_t , Function > > enrichfunc ;
	Form * behaviour ;
	NonLinearForm * nonlinbehaviour ;
	ElementState * state ;
public:
	
	ElementarySurface(bool f = false) ;
	virtual Function jacobian() const = 0;
	virtual ~ElementarySurface() ;
	
	virtual void print()  const = 0 ;
	
	const bool isFather ;
	
	virtual bool isMoved() const = 0 ;
	
	virtual const std::valarray< Function > & getShapeFunctions() const ;

	virtual const std::vector<std::pair<size_t, Function> > getDofs() const ;
	virtual const std::vector< size_t > getDofIds() const ;
	
	virtual const Function & getShapeFunction(size_t i) const ;
	virtual Function & getShapeFunction(size_t i)  ;
	
	virtual const std::pair<size_t, Function> & getEnrichmentFunction(size_t i) const ;
	virtual std::pair<size_t, Function> & getEnrichmentFunction(size_t i)  ;
	
	virtual const std::vector< std::pair< size_t, Function> > & getEnrichmentFunctions() const ;
	
	virtual  std::vector< std::pair< size_t, Function> >  & getEnrichmentFunctions()  ;
	virtual const Function getXTransform() const ;
	virtual const Function getYTransform() const ;
	virtual const Function getTTransform() const ;
	const Function getdXTransform(Variable) const ;
	const Function getdYTransform(Variable) const ;
	const Function getdTTransform(Variable) const ;
	const double getdXTransform(Variable, const Point p) const ;
	const double getdYTransform(Variable, const Point p) const ;
	const double getdTTransform(Variable, const Point p) const ;
	
	virtual ElementState * getState() const ;
	
	void setEnrichment(std::pair<size_t,Function> p) ;
	virtual Point inLocalCoordinates(const Point & p) const  = 0;
	
	virtual Matrix getInverseJacobianMatrix(const Point & p) const = 0 ;
	
	virtual std::vector<std::vector<Matrix> > getElementaryMatrix() const = 0;
	virtual Vector getForces() const { return Vector(0) ;}
	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix() const = 0;
	virtual Vector getNonLinearForces() const = 0 ;
	
	virtual const Point &  getPoint(size_t i) const = 0 ;
	virtual  Point &  getPoint(size_t i)  = 0 ;
	
	virtual Form * getBehaviour() const ;
	void setBehaviour(Form *);
	
	virtual NonLinearForm * getNonLinearBehaviour() const;
	void setNonLinearBehaviour(NonLinearForm * f) ;
	
	void step(double dt, Vector * displacements) ;
	void nonLinearStep(double dt, Vector *displacements) ;
	
	Order getOrder() const;
	void setOrder(Order) ;
	
} ;


class TriElement : public Triangle, public ElementarySurface
{
	
protected :
	
	std::valarray< std::pair<Point, double> > genGaussPoints() const;
	virtual void computeCenter();
	
public:
	
	bool moved ;
	
	GEO_DERIVED_OBJECT(Triangle) ;
	
	TriElement( Point * p0,  Point * p1,  Point * p2, bool father = false) ;
	
	TriElement(Order order = LINEAR, bool father = true) ;
	void refresh(const TriElement * parent) ;
	
	virtual std::vector<std::vector<Matrix> > getElementaryMatrix() const ;
	
	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix(Vector * state) const ;
	
	Function jacobian() const ;
	
	double  jacobianAtPoint(const Point p) const ;
	
	Matrix getInverseJacobianMatrix(const Point & p) const ;
	
	std::valarray< std::pair<Point, double> > getGaussPoints() const;
	
	virtual bool isMoved() const;
	
	virtual void print() const;
	
	virtual Point inLocalCoordinates(const Point &p) const ;
	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix() const ;
	
	virtual Vector getNonLinearForces() const ;
	
} ;


class ElementaryVolume : public IntegrableEntity
{
protected:
	std::valarray< Function > *shapefunc ;
	std::vector< std::pair<size_t , Function> > enrichfunc ;
	virtual std::valarray< std::pair<Point, double> > genGaussPoints() const = 0 ;
	Form * behaviour ;
	Order order ;
	ElementState * state ;
	NonLinearForm * nonlinbehaviour ;
	
public:
	ElementaryVolume(bool f = false) ;
	virtual Function jacobian() const ;
	double jacobianAtPoint(const Point & p) const ;
	virtual ~ElementaryVolume() ;
	
	const bool isFather ;
	
	virtual bool isMoved() const ;
	
	virtual const std::vector<std::pair<size_t, Function> > getDofs() const ;
	virtual const std::vector< size_t > getDofIds() const ;
	
	virtual void print()  const = 0 ;
	virtual std::valarray< std::pair<Point, double> > getGaussPoints() const { return genGaussPoints() ;}
	virtual std::vector<std::vector<Matrix> > getElementaryMatrix() const = 0;
	virtual Form * getBehaviour() const ;
	virtual void setBehaviour(Form *);
	virtual Vector getForces() const { return Vector(0) ;}

	virtual NonLinearForm * getNonLinearBehaviour() const ;
	
	virtual const std::valarray< Function > & getShapeFunctions() const ;
	virtual const Function & getShapeFunction(size_t i) const ;
	virtual Function & getShapeFunction(size_t i)  ;
	virtual const std::pair<size_t, Function> & getEnrichmentFunction(size_t i) const  ;
	virtual  std::pair<size_t, Function> & getEnrichmentFunction(size_t i)   ;
	virtual const std::vector< std::pair< size_t, Function> >  & getEnrichmentFunctions() const ;
	virtual std::vector< std::pair< size_t, Function> >  & getEnrichmentFunctions()  ;
	virtual const Function getXTransform() const ;
	virtual const Function getYTransform() const ;
	virtual const Function getZTransform() const ;
	virtual const Function getTTransform() const ;
	virtual const Function getdXTransform(Variable v) const ;
	virtual const Function getdYTransform(Variable v) const ;
	virtual const Function getdZTransform(Variable v) const ;
	virtual const Function getdTTransform(Variable v) const ;
	virtual const double getdXTransform(Variable v, const Point & p) const ;
	virtual const double getdYTransform(Variable v, const Point & p) const ;
	virtual const double getdZTransform(Variable v, const Point & p) const ;
	virtual const double getdTTransform(Variable v, const Point & p) const ;

	virtual void setEnrichment(std::pair<size_t, Function>  p) ;
	
	virtual Matrix getInverseJacobianMatrix(const Point & p) const ;
	
	virtual ElementState * getState() const ;
	
	virtual const std::valarray< Point * > & getBoundingPoints() const = 0;
	virtual std::valarray< Point * > & getBoundingPoints() = 0;
	virtual const Point & getBoundingPoint(size_t i) const = 0;
	virtual Point & getBoundingPoint(size_t i) = 0;
	virtual const std::valarray< Point * > & getInPoints() const = 0;
	virtual std::valarray< Point * > & getInPoints() = 0;
	virtual const Point &  getPoint(size_t i) const = 0 ;
	virtual Point &  getPoint(size_t i) = 0 ;
		
	
	virtual void step(double dt, Vector * displacements) ;
	virtual void nonLinearStep(double dt, Vector *displacements) ;
	
	virtual Order getOrder() const;
	
} ;


class TetrahedralElement : public Tetrahedron,  public ElementaryVolume
{
protected :
	std::valarray< std::pair<Point, double> > genGaussPoints() const;
	virtual void computeCenter();
public:
	bool moved;
	GEO_DERIVED_OBJECT(Tetrahedron) ;
	
	TetrahedralElement( Point * p0,  Point * p1,  Point * p2, Point * p3, bool father = false) ;
	TetrahedralElement( Point * p0,  Point * p1,  Point * p2, Point * p3, Point * p4,  Point * p5,  Point * p6, Point * p7, bool father = false) ;
	TetrahedralElement(Order order = LINEAR, bool father = true);
	TetrahedralElement(TetrahedralElement * parent, Tetrahedron * t);
	virtual std::vector<std::vector<Matrix> > getElementaryMatrix() const ;
	
	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix() const ;
	
	virtual Vector getForces() const ;
	
	virtual Vector getNonLinearForces() const ;
	
	void refresh(const TetrahedralElement * parent);

	virtual void print() const;
	virtual Point inLocalCoordinates(const Point & p) const ;

} ;

class HexahedralElement : public Hexahedron,  public ElementaryVolume
{
protected :
	std::valarray< std::pair<Point, double> > genGaussPoints() const;
	virtual void computeCenter();
public:
	std::vector<HexahedralElement *> neighbourhood ;
// 	bool moved;
	GEO_DERIVED_OBJECT(Hexahedron) ;
	
	HexahedralElement(Order order, bool f = true) ;
	HexahedralElement(HexahedralElement * parent,Hexahedron * t);

	virtual std::vector<std::vector<Matrix> > getElementaryMatrix() const ;

	virtual Vector getForces() const ;

	void refresh(const HexahedralElement * parent);
	virtual void print()  const;
	virtual Point inLocalCoordinates(const Point & p) const ;
	virtual Vector getNonLinearForces() const;
	virtual std::vector<std::vector<Mu::Matrix> > getNonLinearElementaryMatrix() const;
	bool visited ;
} ;

} ;

void computeNeighbourhoodForStructuredHexahedralMesh(std::vector<Mu::HexahedralElement *> & vec) ;
void burn(std::vector<Mu::HexahedralElement *> & vec) ;

#endif // __ELEMENTS_H_
