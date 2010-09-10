
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef INTERGRABLE_ENTITY_H
#define INTERGRABLE_ENTITY_H

#include "../geometry/geometry_base.h"
#include "../polynomial/vm_base.h"
#include "../polynomial/vm_function_matrix.h"
#include "../utilities/xml.h"
#include "../physics/homogenization/homogenization_base.h"

namespace Mu
{

class DelaunayTreeItem3D;
class DelaunayTetrahedron;
class DelaunayTreeItem;
class Assembly ;

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
	QUADTREE_REFINED,
	REGULAR_GRID
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
struct FractureCriterion ;
struct VirtualMachine ;
struct BoundaryCondition ;

/** \brief State of the element, allows easy extraction of the various fields
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
	
	/** \brief Construct the state of the argument*/
	ElementState(IntegrableEntity *) ;
	/** \brief Copy-constructor*/
	ElementState(const ElementState &s) ;
						
	ElementState & operator =(const ElementState &) ;
	
/** \brief Return strain at given point*/
	Vector getStrain(const Point & , bool local = false) const;

/** \brief Return stress at given point*/
	Vector getPreviousStress(const Point & , bool local = false) const;
	
/** \brief Return stress at given point*/
	Vector getStress(const Point & , bool local = false) const;

/** \brief Return stress at given point, ignoring enrichment functions*/
	Vector getNonEnrichedStress(const Point & , bool local = false) const;

/** \brief Return strain at given point, ignoring enrichment functions*/
	Vector getNonEnrichedStrain(const Point & , bool local = false) const;

/** \brief Return strain at given point, in matrix form*/
	Matrix getStrainMatrix(const Point & , bool local = false) const;

/** \brief Return stress at given point, in matrix form*/
	Matrix getPreviousStressMatrix(const Point & p, bool local = false) const;
	
/** \brief Return stress at given point, in matrix form*/
	Matrix getStressMatrix(const Point & , bool local = false) const;

/** \brief Return strain at given point*/
	Vector getStrain(const std::pair<Point, double> &  ) const;

/** \brief Return stress at given point*/
	Vector getStress(const std::pair<Point, double> &  ) const;

/** \brief Return strain at given points*/
	Vector getStrain(const Mu::PointArray &) const;

/** \brief Return stress at given points*/
	Vector getPreviousStress(const Mu::PointArray & v) const;
	
/** \brief Return stress at given points*/
	Vector getStress(const Mu::PointArray &) const;

/** \brief Return stress at given points, ignoring enrichment functions*/
	Vector getNonEnrichedStress(const Mu::PointArray &) const;

/** \brief Return strain at given points, ignoring enrichment functions*/
	Vector getNonEnrichedStrain(const Mu::PointArray &) const;

/** \brief Return strain at given (Gauss) points*/
	Vector getStrain(const std::valarray<std::pair<Point, double> > &p) const;

/** \brief Return stress at given (Gauss) points*/
	Vector getStress(const std::valarray<std::pair<Point, double> > &p) const;

/** \brief Return stress at given (Gauss) points, ignoring enrichment functions*/
	Vector getNonEnrichedStress(const std::valarray<std::pair<Point, double> > & p) const;

/** \brief Return stress at given (Gauss) points, ignoring enrichment functions, given inverse jacobian matrices*/
	Vector getNonEnrichedStress(const std::valarray<std::pair<Point, double> > & p, const std::valarray<Matrix> &Jinv) const;

/** \brief Return strain at given (Gauss) points, ignoring enrichment functions*/
	Vector getNonEnrichedStrain(const std::valarray<std::pair<Point, double> > & p) const;

/** \brief Return strain at given (Gauss) points, ignoring enrichment functions, given inverse jacobian matrices*/
	Vector getNonEnrichedStrain(const std::valarray<std::pair<Point, double> > & p, const std::valarray<Matrix> &Jinv) const;

/** \brief return stress and strain at given points*/
	std::pair<Vector, Vector > getStressAndStrain(const Mu::PointArray &) const;

/** \brief return stress and strain at given (gauss) points*/
	std::pair<Vector, Vector > getStressAndStrain( std::valarray<std::pair<Point, double> > & p) const;
	
/** \brief get Principal Stresses at given point*/
	Vector getPrincipalStresses(const Point & , bool local = false) const ;
	
/** \brief get Principal Stresses at given point*/
	Vector getPreviousPrincipalStresses(const Point & , bool local = false) const ;

/** \brief get Principal Stresses at given points*/
	Vector getPrincipalStresses(const Mu::PointArray &) const ;
	
/** \brief get Principal Stresses at given points*/
	Vector getPreviousPrincipalStresses(const Mu::PointArray &) const ;

/** \brief get Principal Strains at given point*/
	Vector getPrincipalStrains(const Point & p, bool local = false) const ;
	
/** \brief get Principal Strains at given points*/
	Vector getPrincipalStrains(const Mu::PointArray &) const ;
	
/** \brief get Principal angle at given point*/
	double getPrincipalAngle(const Point & p, bool local= false) const ;

/** \brief get Principal angle at given points*/
	Vector getPrincipalAngle(const Mu::PointArray & v) const;
	
/** \brief get symbolic expression of Stress, given he inverse Jacobian*/
	FunctionMatrix getStressFunction(const Matrix &Jinv) const;

/** \brief get symbolic expression of Strain, given he inverse Jacobian*/
	FunctionMatrix getStrainFunction(const Matrix &Jinv) const;

/** \brief get symbolic expression of displacement, given he inverse Jacobian*/
	FunctionMatrix getDisplacementFunction() const;
	
/** \brief get average desplacements over the element*/
	Vector getAverageDisplacement() const ;
	
	double getVonMisesStrain(const Mu::Point& p, bool local = false) const ;
	
/** \brief return maximum Von Mises Stress*/
	double getMaximumVonMisesStress() const ;
	
/** \brief return maximum Von Mises Stress*/
	double getPreviousMaximumVonMisesStress() const ;
	
/** \brief return displacement at point*/
	Vector getDisplacements(const Point &, bool local = false, bool fast = false, const Vector * source = NULL) const;

/** \brief return the linear interpolating factors for the displacement field at the given point*/
	std::vector<double> getInterpolatingFactors(const Point & p, bool local = false) const ;

/** \brief return the linear enrichment interpolating factors for the displacement field at the given point*/
	std::vector<double> getEnrichedInterpolatingFactors(const Point & p, bool local = false) const ;
	
/** \brief return the linear non-enrichment interpolating factors for the displacement field at the given point*/
	std::vector<double> getNonEnrichedInterpolatingFactors(const Point & p, bool local = false) const ;
	
/** \brief return displacement at points*/
	Vector getDisplacements(const std::valarray<Point> & p) const ;

/** \brief return displacement at points*/
	Vector getDisplacements(const std::vector<std::pair<Point, double> > & p ) const ;

/** \brief return displacements at the nodes of the element*/
	const Vector & getDisplacements() const;

/** \brief Return the elastic energy of this element*/
	double elasticEnergy() const ;

/** \brief return displacements at the nodes of the element*/
	Vector & getDisplacements() ;
	
/** \brief return previous displacements at point*/
	Vector getPreviousDisplacements(const Point &) const;

/** \brief return previous displacements at points*/
	Vector getPreviousDisplacements(const std::valarray<Point> & p) const ;

/** \brief return previous displacements at the nodes of the element*/
	const Vector & getPreviousDisplacements() const;

/** \brief return previous displacements at the nodes of the element*/
	Vector & getPreviousDisplacements() ;
	
/** \brief return displacements two timesteps back at point*/
	Vector  getPreviousPreviousDisplacements(const Point &) const;

/** \brief return displacements two timesteps back at point*/
	Vector getPreviousPreviousDisplacements(const std::valarray<Point> & p) const ;

/** \brief return displacements two timesteps back at the nodes of the element*/
	const Vector & getPreviousPreviousDisplacements() const;

/** \brief return displacements two timesteps back at the nodes of the element*/
	Vector & getPreviousPreviousDisplacements() ;
	
/** \brief get enriched displacements at the nodes of the element*/
	const Vector & getEnrichedDisplacements() const;

/** \brief get enriched displacements at the nodes of the element*/
	Vector & getEnrichedDisplacements() ;

/** \brief get previous enriched displacements at the nodes of the element*/
	Vector & getPreviousEnrichedDisplacements() ;

/** \brief get previous enriched displacements at the nodes of the element*/
	const Vector & getPreviousEnrichedDisplacements() const;

/** \brief get enriched displacements at the nodes of the element from two timesteps back*/
	const Vector & getPreviousPreviousEnrichedDisplacements() const;

/** \brief get enriched displacements at the nodes of the element from two timesteps back*/
	Vector & getPreviousPreviousEnrichedDisplacements() ;
	
	void stepBack() ;
	
	Vector getDeltaStrain(const Point & ) const;
	Vector getDeltaStress(const Point & ) const;
	Vector getDeltaStrain(const std::pair<Point, double> &  ) const;
	Vector getDeltaStress(const std::pair<Point, double> &  ) const;
	Vector getDeltaStrain(const Mu::PointArray &) const;
	Vector getDeltaStress(const Mu::PointArray &) const;
	Vector getDeltaStrain(const std::valarray<std::pair<Point, double> > & p) const;
	Vector getDeltaStress(const std::valarray<std::pair<Point, double> > & p) const;
	std::pair<Vector, Vector > getDeltaStressAndDeltaStrain(const Mu::PointArray &) const;
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
	
	void step(double dt, const Vector* d ) ;
	
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

/** \brief container for a set of Gauss points*/
struct GaussPointArray
{
	std::valarray< std::pair<Point, double> > gaussPoints ;
	int id ;
	GaussPointArray() : gaussPoints(std::make_pair(Point(), 1.),1), id(-2) { } ;
	GaussPointArray(const std::valarray< std::pair<Point, double> > & array, int i): gaussPoints(array), id(i) { } ;
	void operator = (const GaussPointArray & gp) 
	{
		gaussPoints.resize(gp.gaussPoints.size()) ; 
		gaussPoints = gp.gaussPoints ; 
		id = gp.id ;
	}
} ;

/** \brief Abstract class for the representation of elements
 */
class IntegrableEntity : public Geometry
{
protected:
	
	Order order ;
	ElementState state ;
	std::vector<BoundaryCondition *> * boundaryConditionCache ;
	GaussPointArray * cachedGps ;
	const GaussPointArray * getCachedGaussPoints() const {return cachedGps ;};
	void setCachedGaussPoints(GaussPointArray * gp) { cachedGps = gp ;};
public:
	bool enrichmentUpdated ;
	bool behaviourUpdated ;
	bool enrichmentFunctionsCompiled ;
	
	IntegrableEntity() ;
	virtual void getInverseJacobianMatrix(const Point &p, Matrix & ret) const = 0 ;
	virtual ~IntegrableEntity() ;
	virtual const Point & getPoint(size_t i) const = 0 ;
	virtual Point & getPoint(size_t i)  = 0 ;
	virtual const GaussPointArray & getGaussPoints() = 0;
	virtual Point inLocalCoordinates(const Point & p) const  = 0;
	virtual double area() const { return 0 ; } 
	virtual double volume() const { return 0 ; } 
	
	virtual Function getXTransform() const = 0;
	virtual Function getYTransform() const = 0;

	virtual	const std::valarray< Function >  & getShapeFunctions() const = 0 ;
	virtual	const std::vector< Function> & getEnrichmentFunctions() const = 0 ;
// 	virtual	std::vector< Function> & getEnrichmentFunctions() = 0 ;
	virtual	const Function & getShapeFunction(size_t i) const = 0 ;
// 	virtual	Function & getShapeFunction(size_t i)  = 0 ;
	virtual const Function & getEnrichmentFunction(size_t i) const = 0;	
// 	virtual Function & getEnrichmentFunction(size_t i) = 0;	
	virtual Order getOrder() const  = 0 ;
	
	virtual void compileAndPrecalculate() = 0 ;
	virtual std::vector<size_t> clearEnrichment(const Geometry * g) = 0 ;
	virtual const std::vector< size_t > getDofIds() const = 0;
	
	virtual Form * getBehaviour() const = 0;
	virtual NonLinearForm * getNonLinearBehaviour() const = 0;
	virtual std::vector<std::vector<Matrix> > & getElementaryMatrix()  = 0 ;
	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix()  = 0 ;
// 	virtual Vector getForces() = 0 ;
	virtual Vector getNonLinearForces() = 0 ;
	virtual void applyBoundaryCondition(Assembly * a) ;
	
	virtual bool isMoved() const = 0 ;
	virtual void print() const = 0;
	
	virtual const ElementState & getState() const ;
	virtual ElementState & getState() ;
	
	virtual Mesh<DelaunayTriangle, DelaunayTreeItem> * get2DMesh() const = 0 ;
	virtual Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * get3DMesh() const = 0 ;
	
} ;

/** \brief A Form for DIM degrees of freedom.
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
	 * @param e The element in which the integration is performed
	 * @return The vector of the symbolic matrices at the integration points
	 */
	virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &, VirtualMachine * vm) const = 0 ;

	virtual XMLTree * toXML() {return new XMLTree("abstract form") ; } ;
	virtual Material toMaterial()
	{
		Material mat(getTensor(Point(0,0))) ;
		return mat ;
	}
	
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
	
	virtual bool hasInducedBoundaryConditions() const
	{
		return false ;
	} ;
	
	virtual size_t getNumberOfDegreesOfFreedom() const { return this->num_dof ; }
	
	virtual void transform(const Function & x, const Function & y) { } ;
	virtual void transform(const Function & x, const Function & y, const Function & z) { } ;
	
	/** Step through time
	 * 
	 * @param timestep length of the timestep
	 * @param currentState State of the element with this behaviour
	 */
	virtual void step(double timestep, ElementState & currentState) = 0;
	virtual void artificialDamageStep(double d) = 0 ;
	virtual void updateElementState(double timestep, ElementState & currentState) const = 0;
	
	virtual bool fractured() const = 0 ;
	virtual bool changed() const { return false ; } ;
	virtual void getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & v) const { };
	virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const { 
		return  std::vector<BoundaryCondition * >() ;
	} ;
	
	virtual Form * getCopy() const = 0 ;
	virtual void stepBack() { }  ;
	
	virtual Matrix getTensor(const Point & p) const
	{
		return param ;
	}
	
	virtual Matrix getPreviousTensor(const Point & p) const
	{
		return param ;
	}

	virtual Vector getImposedStress(const Point & p) const
	{
		return Vector(double(0), param.numCols()) ;
	}
	
	virtual ~Form() { } ;

	virtual FractureCriterion * getFractureCriterion() const { return NULL ; }
} ;

} ;


#endif
