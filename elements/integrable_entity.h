
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
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
class FractureCriterion ;

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

typedef enum{
	REAL_STRESS,
	EFFECTIVE_STRESS
} StressCalculationMethod ; 

typedef enum {
	DISPLACEMENT_FIELD,
	ENRICHED_DISPLACEMENT_FIELD,
	SPEED_FIELD,
	FLUX_FIELD,
	GRADIENT_FIELD,
	STRAIN_FIELD,
	STRAIN_RATE_FIELD,
	EFFECTIVE_STRESS_FIELD,
	REAL_STRESS_FIELD,
	PRINCIPAL_STRAIN_FIELD,
	PRINCIPAL_EFFECTIVE_STRESS_FIELD,
	PRINCIPAL_REAL_STRESS_FIELD,
	NON_ENRICHED_STRAIN_FIELD,
	NON_ENRICHED_STRAIN_RATE_FIELD,
	NON_ENRICHED_EFFECTIVE_STRESS_FIELD,
	NON_ENRICHED_REAL_STRESS_FIELD,
	VON_MISES_STRAIN_FIELD,
	VON_MISES_REAL_STRESS_FIELD,
	VON_MISES_EFFECTIVE_STRESS_FIELD,
	PRINCIPAL_ANGLE_FIELD,
	INTERNAL_VARIABLE_FIELD,
	GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD,
	GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD,
	GENERALIZED_VISCOELASTIC_SPEED_FIELD,
	GENERALIZED_VISCOELASTIC_STRAIN_FIELD,
	GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD,
	GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD,
	GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD,
	GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD,
	GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD,
	GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD,
	GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD,
	GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD,
	GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD,
	GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD,	
} FieldType ;

struct Form ;
struct NonLinearForm ;
struct Function ;
struct DelaunayTriangle ;
struct IntegrableEntity ;
struct FractureCriterion ;
struct DamageModel ;
struct VirtualMachine ;
struct BoundaryCondition ;

/** \brief State of the element, allows easy extraction of the various fields
 * 
 */
class ElementState
{
protected:

	Vector strainAtNodes ;
	Vector stressAtNodes ;
	Vector effectivePStressAtGaussPoints ;

	Vector strainAtGaussPoints ;
	Vector stressAtGaussPoints ;
	
	Vector pstrainAtGaussPoints ;
	Vector pstressAtGaussPoints ;
	
	Vector displacements ;
	Vector enrichedDisplacements ;
	
	Vector previousDisplacements ;
	Vector previousEnrichedDisplacements ;

	
	Vector buffer ;
	
	double timePos ;
	double previousTimePos ;
	
	IntegrableEntity * parent ;
	
	std::vector<ElementState> history ;
	
public:
	
	/** \brief Construct the state of the argument*/
	ElementState(IntegrableEntity *) ;
	/** \brief Copy-constructor*/
	ElementState(const ElementState &s) ;
						
	ElementState & operator =(const ElementState &) ;
	
	virtual void getExternalField( Vector & nodalValues, int externaldofs, const Point & p, Vector & ret, bool local) const ;

	virtual void getExternalFieldAtGaussPoints( Vector & nodalValues, int externaldofs, std::vector<Vector> & ret) const ;
	
	virtual void getField( FieldType f, const Point & p, Vector & ret, bool local, int i = 0) const ;
		
	virtual void getField( FieldType f, const PointArray & p, Vector & ret, bool local, int i = 0) const  ;

	virtual void getField( FieldType f, const std::valarray<std::pair<Point, double> > & p, Vector & ret, bool local, int i = 0) const  ;

	virtual void getFieldAtNodes( FieldType f, Vector & ret, int i = 0) ;

	virtual void getFieldAtGaussPoint( FieldType f, size_t g, Vector & ret, int i = 0) ;
	
	virtual void getField( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, int i = 0, int j = 0) const  ;
		
	virtual void getField( FieldType f1, FieldType f2, const PointArray & p, Vector & ret1, Vector & ret2, bool local, int i = 0, int j = 0) const  ;

	virtual void getField( FieldType f1, FieldType f2, const std::valarray<std::pair<Point, double> > & p, Vector & ret1, Vector & ret2, bool local, int i = 0, int j = 0) const  ;

	virtual void getFieldAtNodes( FieldType f1, FieldType f2, Vector & ret1, Vector & ret2, int i = 0, int j = 0) ;

	virtual void getFieldAtGaussPoint( FieldType f1, FieldType f2, size_t g, Vector & ret1, Vector & ret2, int i = 0, int j = 0) ;
	
	virtual void getAverageField( FieldType f, Vector & ret, int i= 0) ;
	

/** \brief return displacements at the nodes of the element*/
	const Vector & getDisplacements() const;

/** \brief return displacements at the nodes of the element*/
	Vector & getDisplacements() ;

/** \brief return enriched displacements at the nodes of the element*/
	const Vector & getEnrichedDisplacements() const;

/** \brief return enriched displacements at the nodes of the element*/
	Vector & getEnrichedDisplacements() ;

/** \brief return displacements at the nodes of the element*/
	const Vector & getPreviousDisplacements() const;

/** \brief return displacements at the nodes of the element*/
	Vector & getPreviousDisplacements() ;

/** \brief return enriched displacements at the nodes of the element*/
	const Vector & getPreviousEnrichedDisplacements() const;

/** \brief return enriched displacements at the nodes of the element*/
	Vector & getPreviousEnrichedDisplacements() ;

		
/** \brief Return the set of eigenvector forming the reference frame of the principal stresses*/
	std::vector<Point> getPrincipalFrame( const Point &p, bool local = false, StressCalculationMethod m = REAL_STRESS)  ;
	
/** \brief get symbolic expression of Stress, given he inverse Jacobian*/
	FunctionMatrix getStressFunction(const Matrix &Jinv, StressCalculationMethod m = REAL_STRESS) const;

/** \brief get symbolic expression of Strain, given he inverse Jacobian*/
	FunctionMatrix getStrainFunction(const Matrix &Jinv) const;

/** \brief get symbolic expression of displacement, given he inverse Jacobian*/
	FunctionMatrix getDisplacementFunction() const;
	
/** \brief get average desplacements over the element*/
	Vector getAverageDisplacement() const ;
	
/** \brief return the linear interpolating factors for the displacement field at the given point*/
	std::vector<double> getInterpolatingFactors(const Point & p, bool local = false) const ;

/** \brief return the linear enrichment interpolating factors for the displacement field at the given point*/
	std::vector<double> getEnrichedInterpolatingFactors(const Point & p, bool local = false) const ;
	
/** \brief return the linear non-enrichment interpolating factors for the displacement field at the given point*/
	std::vector<double> getNonEnrichedInterpolatingFactors(const Point & p, bool local = false) const ;
	
/** \brief Return the elastic energy of this element*/
	double elasticEnergy() ;

	virtual void step(double dt, const Vector* d ) ;
	virtual void elasticStep(double dt, const Vector* d ) { } ;
	
	double getTime() const ;
	double getDeltaTime() const ;
	
	IntegrableEntity * getParent() const
	{
		return parent ;
	}
	
	virtual void initialize(bool initializeFractureCache =true) ;
	
	const Vector & getBuffer() const ;
	Vector & getBuffer()  ;
	
} ;

class ElementStateWithInternalVariables : public ElementState
{
protected:
	std::vector<std::vector<Vector> > internalVariablesAtGaussPoints ;
	int n ;
	int p ;
	
public:
	ElementStateWithInternalVariables(IntegrableEntity *, int n, int p) ;

	ElementStateWithInternalVariables(const ElementStateWithInternalVariables &s) ;
						
	ElementStateWithInternalVariables & operator =(const ElementStateWithInternalVariables &) ;
	
	int numberOfInternalVariables() const { return n ; }
	
	int sizeOfInternalVariable() const { return p ; }
	
	virtual void getField( FieldType f, const Point & p, Vector & ret, bool local, int i = 0) const ;
		
	virtual void getField( FieldType f, const PointArray & p, Vector & ret, bool local, int i = 0) const  ;

	virtual void getField( FieldType f, const std::valarray<std::pair<Point, double> > & p, Vector & ret, bool local, int i = 0) const  ;
		
	virtual void getFieldAtNodes( FieldType f, Vector & ret, int i = 0) ;

	virtual void getFieldAtGaussPoint( FieldType f, size_t g, Vector & ret, int i = 0) ;

	virtual void getField( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, int i = 0, int j = 0) const  ;
		
	virtual void getField( FieldType f1, FieldType f2, const PointArray & p, Vector & ret1, Vector & ret2, bool local, int i = 0, int j = 0) const  ;

	virtual void getField( FieldType f1, FieldType f2, const std::valarray<std::pair<Point, double> > & p, Vector & ret1, Vector & ret2, bool local, int i = 0, int j = 0) const  ;

	virtual void getFieldAtNodes( FieldType f1, FieldType f2, Vector & ret1, Vector & ret2, int i = 0, int j = 0) ;

	virtual void getFieldAtGaussPoint( FieldType f1, FieldType f2, size_t g, Vector & ret1, Vector & ret2, int i = 0, int j = 0) ;

	virtual void initialize(bool initializeFractureCache =true) ;
	
	virtual void setInternalVariableAtGaussPoint(Vector & v, size_t g, int i) ;
	
} ;

class ParallelElementState : public ElementState
{
protected:
	std::vector<ElementState *> states ;

public:	
	ParallelElementState(IntegrableEntity *, std::vector<ElementState *> s) ;
	ParallelElementState(const ParallelElementState &s) ;						
	ParallelElementState & operator =(const ParallelElementState & s) ;
	
	ElementState & getState(size_t i) ;
	const ElementState & getState(size_t i) const ;
	size_t getNumberOfStates() const { return states.size() ; }

	virtual void initialize(bool initializeFractureCache =true) ;
	virtual void step(double dt, const Vector* d ) ;
	
} ;

class SerialElementState : public ElementState
{
protected:
	std::vector<ElementState *> states ;

public:	
	SerialElementState(IntegrableEntity *, std::vector<ElementState *> s) ;
	SerialElementState(const SerialElementState &s) ;						
	SerialElementState & operator =(const SerialElementState & s) ;
	
	ElementState & getState(size_t i) ;
	const ElementState & getState(size_t i) const ;
	size_t getNumberOfStates() const { return states.size() ; }

	virtual void initialize(bool initializeFractureCache =true) ;
	virtual void step(double dt, const Vector* d ) ;
	
} ;

class KelvinVoightSpaceTimeElementState : public ElementState
{
public:
	KelvinVoightSpaceTimeElementState(IntegrableEntity * e) ;
	KelvinVoightSpaceTimeElementState(const KelvinVoightSpaceTimeElementState &s) ;						
	KelvinVoightSpaceTimeElementState & operator =(const KelvinVoightSpaceTimeElementState & s) ;
	
	virtual void getField( FieldType f, const Point & p, Vector & ret, bool local, int i = 0) const ;
		
	virtual void getField( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, int i = 0, int j = 0) const  ;
			
	virtual void getFieldAtNodes( FieldType f1, FieldType f2, Vector & ret1, Vector & ret2, int i = 0, int j = 0) ;

} ;


/** \brief container for a set of Gauss points*/
struct GaussPointArray
{
	std::valarray< std::pair<Point, double> > gaussPoints ;
	int id ;
	GaussPointArray() : gaussPoints(std::make_pair(Point(), 1.),1), id(-2) { } ;
	GaussPointArray(const std::pair<Point, double> & p) : gaussPoints(p, 1), id(-2) { } ;
	GaussPointArray(const GaussPointArray & gp) : gaussPoints(gp.gaussPoints), id(-2) { } ;
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
	ElementState * state ;
	std::vector<BoundaryCondition *> * boundaryConditionCache ;
	GaussPointArray * cachedGps ;
	const GaussPointArray * getCachedGaussPoints() const { return cachedGps ;};
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
	virtual Function getZTransform() const ;

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
	virtual std::vector<size_t> clearAllEnrichment() = 0;
	virtual const std::vector< size_t > getDofIds() const = 0;
	
	virtual Form * getBehaviour() const = 0;
	virtual NonLinearForm * getNonLinearBehaviour() const = 0;
	virtual std::valarray<std::valarray<Matrix> > & getElementaryMatrix()  = 0 ;
	virtual std::valarray<std::valarray<Matrix> > getNonLinearElementaryMatrix()  = 0 ;
// 	virtual Vector getForces() = 0 ;
	virtual Vector getNonLinearForces() = 0 ;
	virtual void applyBoundaryCondition(Assembly * a) ;
	
	virtual bool isMoved() const = 0 ;
	virtual void print() const = 0;
	
	virtual const ElementState & getState() const ;
	virtual ElementState & getState() ;
	virtual void setState( ElementState * state) ;
	
	virtual Mesh<DelaunayTriangle, DelaunayTreeItem> * get2DMesh() const = 0 ;
	virtual Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * get3DMesh() const = 0 ;
	virtual void setOrder( Order o) = 0;
	
} ;

/** \brief A Form for DIM degrees of freedom.
 */
class Form
{
protected:
	bool time_d ;
	bool space_d ;
	size_t num_dof ;
	
	Geometry * source ;
	
public:
	/** A form has at least a parameter, which takes the shape of a Matrix*/
	Matrix param ;
	/** The type helps to know the available parameters and methods of the subclasses*/
	ParametersType type;
	
	Form(const Matrix & p, bool t = false, bool s = false, size_t numdof = 2) : time_d(t), space_d(s), num_dof(numdof), param(p), source(nullptr) { } ;
	Form() : time_d(false), space_d(false), num_dof(0), param(Matrix(0,0)), source(nullptr){ } ;
	
	/** apply the form on a pair of functions
	 * 
	 * @param p_i First form function
	 * @param p_j Second Form Function
	 * @param e The element in which the integration is performed
	 * @return The vector of the symbolic matrices at the integration points
	 */
	virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &, VirtualMachine * vm) const = 0 ;

	virtual XMLTree * toXML() {return new XMLTree("abstract form") ; } ;
	
	virtual bool timeDependent() const
	{
		return this->time_d ;
	} ;
	
	virtual bool spaceDependent() const
	{
		return this->space_d ;
	} ;
	
	virtual Vector getForcesFromAppliedStress( Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v) ;

	virtual Vector getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v) ;
	
	virtual Geometry * getSource() const { return source  ;}
	virtual void setSource( Geometry * const src ) {source = src ;}
	 
	virtual bool hasInducedForces() const ;
	
	virtual void scale(double d) ;
	
	virtual size_t getNumberOfDegreesOfFreedom() const { return num_dof ; }
	
	virtual void transform(const ElementarySurface *) { } ;
	virtual void transform(const ElementaryVolume *) { } ;
	
	/** Step through time
	 * 
	 * @param timestep length of the timestep
	 * @param currentState State of the element with this behaviour
	 */
	virtual void step(double timestep, ElementState & currentState, double maxScore) = 0;
	virtual void updateElementState(double timestep, ElementState & currentState) const = 0;
	
	virtual bool fractured() const = 0 ;
	virtual bool changed() const { return false ; } ;
	virtual void getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & v) const { };
	virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
	
	virtual Form * getCopy() const = 0 ;
	
	virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const
	{
		return param ;
	}
	
	virtual void setTensor(const Matrix & m)
	{
	 param = m ;
	}
	
	virtual Matrix getPreviousTensor(const Point & p) const
	{
		return param ;
	}

	virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
	virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
	
	virtual ~Form() { } ;

	virtual FractureCriterion * getFractureCriterion() const { return nullptr ; }
	
	virtual void setFractureCriterion(FractureCriterion * frac) { }
	
	virtual DamageModel * getDamageModel() const { return nullptr ; }
	
	virtual ElementState * createElementState( IntegrableEntity * e) ;
	
	virtual void preProcess( double timeStep, ElementState & currentState ) { } ;
	
} ;

Matrix makeStressOrStrainMatrix(const Vector & stressOrStrain) ;

int isGaussPoint(const Point & p, IntegrableEntity * e) ;
Vector toPrincipal(const Vector & stressOrStrain) ;
} ;


#endif
