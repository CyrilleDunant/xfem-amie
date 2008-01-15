
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef INTERGRABLE_ENTITY_H
#define INTERGRABLE_ENTITY_H

#include "../geometry/geometry_base.h"
#include "../polynomial/vm_base.h"
#include "../polynomial/vm_function_matrix.h"

namespace Mu
{

typedef enum
{
	CONSTANT,
	LINEAR,
	QUADRATIC,
	CUBIC,
	QUADRIC,
	QUINTIC,
	CONSTANT_TIME_LINEAR,
	CONSTANT_TIME_QUADRATIC,
	LINEAR_TIME_LINEAR,
	LINEAR_TIME_QUADRATIC,
	QUADRATIC_TIME_LINEAR,
	QUADRATIC_TIME_QUADRATIC,
	CUBIC_TIME_LINEAR,
	CUBIC_TIME_QUADRATIC,	
	QUADRIC_TIME_LINEAR,
	QUADRIC_TIME_QUADRATIC,
	QUINTIC_TIME_LINEAR,
	QUINTIC_TIME_QUADRATIC,
} Order ;

typedef enum
{
	PURE_LINEAR,
	LINEAR_AND_CONSTANT,
	NON_LINEAR,
	VOID_BEHAVIOUR
} ParametersType ;

struct Form ;
struct NonLinearForm ;
struct Function ;
struct DelaunayTriangle ;
struct IntegrableEntity ;

/** State of the element, allows easy extraction of the various fields
 * 
 */
class ElementState
{
protected:
	Vector displacements ;
	Vector enrichedDisplacements ;
	
	Vector previousDisplacements ;
	Vector previousEnrichedDisplacements ;
	
	Vector previousPreviousDisplacements ;
	Vector previousPreviousEnrichedDisplacements ;
	
	Vector buffer ;
	
	double timePos ;
	double previousTimePos ;
	double previousPreviousTimePos ;
	
	IntegrableEntity * parent ;
	
	std::vector<ElementState> history ;
	
public:
	
	ElementState(IntegrableEntity *) ;
	ElementState(const ElementState & e) :  displacements(e.getDisplacements()), 
						enrichedDisplacements(e.getEnrichedDisplacements()), 
						previousDisplacements(e.getDisplacements()), 
						previousEnrichedDisplacements(e.getEnrichedDisplacements()), 
						timePos(e.getTime()),
						previousTimePos(e.getTime()-e.getDeltaTime())
	{
		parent = e.getParent() ;
	}
	
	Vector getStrain(const Point & , bool local = false) const;
	Vector getStress(const Point & , bool local = false) const;
	Vector getNonEnrichedStress(const Point & , bool local = false) const;
	Vector getNonEnrichedStrain(const Point & , bool local = false) const;
	Matrix getStrainMatrix(const Point & , bool local = false) const;
	Matrix getStressMatrix(const Point & , bool local = false) const;
	Vector getStrain(const std::pair<Point, double> &  ) const;
	Vector getStress(const std::pair<Point, double> &  ) const;
	Vector getStrain(const std::valarray<Point *> &) const;
	Vector getStress(const std::valarray<Point *> &) const;
	Vector getNonEnrichedStress(const std::valarray<Point *> &) const;
	Vector getNonEnrichedStrain(const std::valarray<Point *> &) const;
	Vector getStrain(const std::valarray<std::pair<Point, double> > &p) const;
	Vector getStress(const std::valarray<std::pair<Point, double> > &p) const;
	Vector getNonEnrichedStress(const std::valarray<std::pair<Point, double> > & p) const;
	Vector getNonEnrichedStress(const std::valarray<std::pair<Point, double> > & p, const std::valarray<Matrix> &Jinv) const;
	Vector getNonEnrichedStrain(const std::valarray<std::pair<Point, double> > & p) const;
	Vector getNonEnrichedStrain(const std::valarray<std::pair<Point, double> > & p, const std::valarray<Matrix> &Jinv) const;
	std::pair<Vector, Vector > getStressAndStrain(const std::valarray<Point *> &) const;
	std::pair<Vector, Vector > getStressAndStrain( std::valarray<std::pair<Point, double> > & p) const;
	
	Vector getPrincipalStresses(const Point & , bool local = false) const ;
	Vector getPrincipalStresses(const std::valarray<Point *> &) const ;
	double getPrincipalAngle(const Point & p, bool local= false) const ;
	Vector getPrincipalAngle(const std::valarray<Point *> & v) const;
	
	FunctionMatrix getStressFunction(const Matrix &Jinv) const;
	FunctionMatrix getStrainFunction(const Matrix &Jinv) const;
	FunctionMatrix getDisplacementFunction() const;
	
	double getMaximumVonMisesStress() const ;
	
	Vector getDisplacements(const Point &, bool local = false) const;
	Vector getDisplacements(const std::valarray<Point> & p) const ;
	Vector getDisplacements(const std::vector<std::pair<Point, double> > & p ) const ;
	const Vector & getDisplacements() const;
	Vector & getDisplacements() ;
	
	Vector getPreviousDisplacements(const Point &) const;
	Vector getPreviousDisplacements(const std::valarray<Point> & p) const ;
	const Vector & getPreviousDisplacements() const;
	Vector & getPreviousDisplacements() ;
	
	Vector  getPreviousPreviousDisplacements(const Point &) const;
	Vector getPreviousPreviousDisplacements(const std::valarray<Point> & p) const ;
	const Vector & getPreviousPreviousDisplacements() const;
	Vector & getPreviousPreviousDisplacements() ;
	
	const Vector & getEnrichedDisplacements() const;
	Vector & getEnrichedDisplacements() ;
	Vector & getPreviousEnrichedDisplacements() ;
	const Vector & getPreviousEnrichedDisplacements() const;
	const Vector & getPreviousPreviousEnrichedDisplacements() const;
	Vector & getPreviousPreviousEnrichedDisplacements() ;
	
	void stepBack() ;
	
	Vector getDeltaStrain(const Point & ) const;
	Vector getDeltaStress(const Point & ) const;
	Vector getDeltaStrain(const std::pair<Point, double> &  ) const;
	Vector getDeltaStress(const std::pair<Point, double> &  ) const;
	Vector getDeltaStrain(const std::valarray<Point *> &) const;
	Vector getDeltaStress(const std::valarray<Point *> &) const;
	Vector getDeltaStrain(const std::valarray<std::pair<Point, double> > & p) const;
	Vector getDeltaStress(const std::valarray<std::pair<Point, double> > & p) const;
	std::pair<Vector, Vector > getDeltaStressAndDeltaStrain(const std::valarray<Point *> &) const;
	std::pair<Vector, Vector > getDeltaStressAndDeltaStrain(const std::valarray<std::pair<Point, double> > & p) const;
	
	Vector getDeltaDisplacements(const Point &) const;
	Vector getDeltaDisplacements(const std::valarray<Point> & p) const ;
	Vector getDeltaDisplacements() const;
	
	Vector getSpeed(const Point &) const ;
	Vector getSpeed(const std::valarray<Point> & p) const ;
	Vector getSpeed() const ;
	
	Vector getAcceleration(const Point &) const ;
	Vector getAcceleration(const std::valarray<Point> & p) const ;
	Vector getAcceleration() const ;
	
	void step(double dt, Vector * ) ;
	
	double getTime() const ;
	double getDeltaTime() const ;
	
	IntegrableEntity * getParent() const
	{
		return parent ;
	}
	
	void initialize() ;
	
	const Vector & getBuffer() const ;
	Vector & getBuffer()  ;
	
} ;


struct GaussPointArray
{
	std::valarray< std::pair<Point, double> > gaussPoints ;
	int id ;
	GaussPointArray() : gaussPoints(std::make_pair(Point(), 1.),1), id(-2) { } ;
	GaussPointArray(const std::valarray< std::pair<Point, double> > & array, int i): gaussPoints(array), id(i) { } ;
} ;

/** Abstract class for the representation of elements
 */
class IntegrableEntity : public Geometry
{
protected:
	
	Order order ;
	ElementState state ;
	
public:
	
	IntegrableEntity() ;
	virtual Matrix getInverseJacobianMatrix(const Point &p) const = 0 ;
	virtual ~IntegrableEntity() { } ;
	virtual const Point & getPoint(size_t i) const = 0 ;
	virtual Point & getPoint(size_t i)  = 0 ;
	virtual GaussPointArray getGaussPoints() const = 0 ;
	virtual Point inLocalCoordinates(const Point & p) const  = 0;
	virtual double area() const { return 0 ; } 
	virtual double volume() const { return 0 ; } 
	
	virtual const Function getXTransform() const = 0;
	virtual const Function getYTransform() const = 0;

	virtual	const std::valarray< Function >  & getShapeFunctions() const = 0 ;
	virtual	const std::vector< Function> & getEnrichmentFunctions() const = 0 ;
	virtual	std::vector< Function> & getEnrichmentFunctions() = 0 ;
	virtual	const Function & getShapeFunction(size_t i) const = 0 ;
	virtual	Function & getShapeFunction(size_t i)  = 0 ;
	virtual const Function & getEnrichmentFunction(size_t i) const = 0;	
	virtual Function & getEnrichmentFunction(size_t i) = 0;	
	virtual Order getOrder() const  = 0 ;
	
// 	virtual const std::vector<std::pair<size_t,const Function &> > getDofs() const = 0;
	virtual const std::vector< size_t > getDofIds() const = 0;
	
	virtual Form * getBehaviour() const = 0;
	virtual NonLinearForm * getNonLinearBehaviour() const = 0;
	virtual std::vector<std::vector<Matrix> > getElementaryMatrix() const  = 0 ;
	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix() const  = 0 ;
	virtual Vector getForces() const = 0 ;
	virtual Vector getNonLinearForces() const = 0 ;
	
	virtual bool isMoved() const = 0 ;
	virtual void print() const = 0;
	
	virtual const ElementState & getState() const ;
	virtual ElementState & getState() ;
	
} ;

/** A Form for DIM degrees of freedom.
 */
class Form
{
protected:
	bool time_d ;
	bool space_d ;
	size_t num_dof ;
public:
	/** A form has at least a parameter, which takes the shape of a Matrix*/
	Matrix param ;
	/** The type helps to know the available parameters and methods of the subclasses*/
	ParametersType type;
	
	Form(const Matrix & p, bool t = false, bool s = false, size_t numdof = 2) : time_d(t), space_d(s), num_dof(numdof), param(p) { } ;
	Form() : time_d(false), space_d(false), num_dof(0), param(Matrix(0,0)){ } ;
	
	/** apply the form on a pair of functions
	 * 
	 * @param p_i First form function
	 * @param p_j Second Form Function
	 * @param Jinv vector of the inverse jacobian matrices at the integration points
	 * @return The vector of the symbolic matrices at the integration points
	 */
	virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const = 0;
	virtual Matrix apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const = 0 ;
	
	virtual bool timeDependent() const
	{
		return this->time_d ;
	} ;
	
	virtual bool spaceDependent() const
	{
		return this->space_d ;
	} ;
	
	virtual bool hasInducedForces() const
	{
		return false ;
	} ;
	
	virtual size_t getNumberOfDegreesOfFreedom() const { return this->num_dof ; }
	
	virtual void transform(const Function & x, const Function & y) { } ;
	virtual void transform(const Function & x, const Function & y, const Function & z) { } ;
	
	/** Step through time
	 * 
	 * @param timestep length of the timestep
	 * @param variables variables affecting the next step.
	 */
	virtual void step(double timestep, ElementState & currentState) = 0;
	virtual void updateElementState(double timestep, ElementState & currentState) const = 0;
	
	virtual bool fractured() const = 0 ;
	virtual bool changed() const { return false ; } ;
	virtual Vector getForces(const ElementState & s, const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const = 0 ;
	
	virtual Form * getCopy() const = 0 ;
	virtual void stepBack() { }  ;
	
	virtual Matrix getTensor(const Point & p) const
	{
		return param ;
	}
	
	virtual ~Form() { } ;
} ;

} ;


#endif
