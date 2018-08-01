
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef INTERGRABLE_ENTITY_H
#define INTERGRABLE_ENTITY_H

#include "../utilities/matrixops.h"
#include "../geometry/geometry_base.h"

#include "../polynomial/vm_base.h"
#include "../polynomial/vm_function_matrix.h"
#include "../utilities/tensor.h"

namespace Amie
{

class LinearForm;

class DelaunayTriangle;
class DelaunayTetrahedron;
class DelaunayTreeItem3D;
class DelaunayTetrahedron;
class DelaunayTreeItem;
class Assembly ;
class FractureCriterion ;

typedef enum : char {
    CONSTANT = 0,
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

typedef enum : char {
    PURE_LINEAR,
    LINEAR_AND_CONSTANT,
    NON_LINEAR,
    VOID_BEHAVIOUR = 3
} ParametersType ;

typedef enum : char {
    REAL_STRESS,
    EFFECTIVE_STRESS
} StressCalculationMethod ;

typedef enum : char {
    DISPLACEMENT_FIELD,
    ENRICHED_DISPLACEMENT_FIELD,
    SPEED_FIELD,
    FLUX_FIELD,
    GRADIENT_FIELD,
    TOTAL_STRAIN_FIELD,
    TOTAL_FINITE_STRAIN_FIELD,
    STRAIN_RATE_FIELD,
    MECHANICAL_STRAIN_FIELD,
    EFFECTIVE_STRESS_FIELD,
    REAL_STRESS_FIELD,
    REAL_FINITE_STRESS_FIELD,
    PRINCIPAL_TOTAL_STRAIN_FIELD,
    PRINCIPAL_MECHANICAL_STRAIN_FIELD,
    PRINCIPAL_EFFECTIVE_STRESS_FIELD,
    PRINCIPAL_REAL_STRESS_FIELD,
    NON_ENRICHED_STRAIN_FIELD,
    NON_ENRICHED_STRAIN_RATE_FIELD,
    NON_ENRICHED_EFFECTIVE_STRESS_FIELD,
    NON_ENRICHED_REAL_STRESS_FIELD,
    VON_MISES_STRAIN_FIELD,
    VON_MISES_REAL_STRESS_FIELD,
    VON_MISES_EFFECTIVE_STRESS_FIELD,
    PRINCIPAL_STRESS_ANGLE_FIELD,
    PRINCIPAL_STRAIN_ANGLE_FIELD,
    INTERNAL_VARIABLE_FIELD,
    IMPOSED_STRESS_FIELD,
    IMPOSED_STRAIN_FIELD,
    PRINCIPAL_IMPOSED_STRAIN_FIELD,
    PRINCIPAL_IMPOSED_STRESS_FIELD,
    SCALAR_DAMAGE_FIELD,
    TENSOR_DAMAGE_FIELD,
    GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD,
    GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD,
    GENERALIZED_VISCOELASTIC_SPEED_FIELD,
    GENERALIZED_VISCOELASTIC_TOTAL_STRAIN_FIELD,
    GENERALIZED_VISCOELASTIC_MECHANICAL_STRAIN_FIELD,
    GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD,
    GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD,
    GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD,
    GENERALIZED_VISCOELASTIC_PRINCIPAL_TOTAL_STRAIN_FIELD,
    GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD,
    GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD,
    GENERALIZED_VISCOELASTIC_NON_ENRICHED_TOTAL_STRAIN_FIELD,
    GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD,
    GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD,
    GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD,
} FieldType ;





class Form ;
class NonLinearForm ;
class PrandtlReussPlasticStrain ;
class PrandtlGrauertPlasticStrain ;
class LinearForm ;
class Function ;
class DelaunayTriangle ;
class ElementarySurface ;
class ElementaryVolume ;
struct IntegrableEntity ;
class FractureCriterion ;
class DamageModel ;
class CollisionDetector ;
class GeometryBasedEffect ;
class VirtualMachine ;
class BoundaryCondition ;
template<class A, class B>
class Mesh  ;
class TwoDCohesiveForces ;



size_t fieldTypeElementarySize ( FieldType f, SpaceDimensionality dim, size_t blocks = 0 ) ;


/** \brief State of the element, allows easy extraction of the various fields
 *
 */
class ElementState
{
  friend TwoDCohesiveForces ;
  friend LinearForm ;
  friend NonLinearForm ;
  friend PrandtlReussPlasticStrain ;
  friend PrandtlGrauertPlasticStrain ;
  friend DelaunayTriangle ;
protected:

    Vector strainAtGaussPoints ;
    Vector stressAtGaussPoints ;

    Vector pstrainAtGaussPoints ;
    Vector pstressAtGaussPoints ;

    Vector displacements ;
    Vector enrichedDisplacements ;
    Vector localExtrapolatedDisplacements ;
    
    bool strainAtGaussPointsSet ;
    bool stressAtGaussPointsSet ;
    bool pstrainAtGaussPointsSet ;
    bool pstressAtGaussPointsSet ;
    
    bool largeDeformations ;

    double timePos ;
    double previousTimePos ;

#ifdef HAVE_OPENMP
    omp_lock_t lock ;
#else
    bool lock ;
#endif
    
    IntegrableEntity * parent ;

    Mesh< DelaunayTriangle, DelaunayTreeItem > * mesh2d ;
    Mesh< DelaunayTetrahedron, DelaunayTreeItem3D > *mesh3d ;
    
   
    
public: 
    void updateInverseJacobianCache(const Point &p) ;
    void getInverseJacobianMatrix(const Point & p, Matrix & ret, bool inital = false) const ;
    
    Matrix * JinvCache ;
    /** \brief Construct the state of the argument*/
    ElementState ( IntegrableEntity * ) ;
    
    virtual ~ElementState() { if(JinvCache) delete JinvCache ;} ;

    virtual void getExternalField ( Vector & nodalValues, int externaldofs, const Point & p, Vector & ret, bool local, VirtualMachine * vm = nullptr )  ;

    virtual void getExternalFieldAtGaussPoints ( Vector & nodalValues, int externaldofs, std::vector<Vector> & ret, VirtualMachine * vm = nullptr )  ;

    virtual void getField ( FieldType f, const Point & p, Vector & ret, bool local, VirtualMachine * vm = nullptr, int i = 0 )  ;

    virtual void getField ( FieldType f, const PointArray & p, Vector & ret, bool local, VirtualMachine * vm = nullptr, int i = 0 )   ;

    virtual void getField ( FieldType f, const std::valarray<std::pair<Point, double> > & p, Vector & ret, bool local, VirtualMachine * vm = nullptr, int i = 0 )   ;

    virtual void getFieldAtGaussPoint ( FieldType f1, FieldType f2, size_t g, Vector & ret1, Vector & ret2, VirtualMachine * vm = nullptr, int i = 0, int j = 0 ) ;

    virtual void getFieldAtGaussPoint ( FieldType f, size_t g, Vector & ret, VirtualMachine * vm = nullptr, int i = 0 ) ;

    virtual void getField ( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm = nullptr, int i = 0, int j = 0 )   ;

    virtual void getField ( FieldType f1, FieldType f2, const PointArray & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm = nullptr, int i = 0, int j = 0 )   ;

    virtual void getField ( FieldType f1, FieldType f2, const std::valarray<std::pair<Point, double> > & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm = nullptr, int i = 0, int j = 0 )   ;

    virtual double getAverageField ( Amie::FieldType f, Vector& ret, Amie::VirtualMachine* vm = nullptr, double t = 0, const std::vector< double > & weights = std::vector<double>(), int index = 0) ;

    virtual double getAverageField ( FieldType f, FieldType f_, Vector& ret, Vector& ret_, VirtualMachine* vm = nullptr, double t = 0, const std::vector< double > & weights = std::vector<double>(), int index = 0 )  ;

    /** \brief return displacements at the nodes of the element*/
    const Vector & getDisplacements() const;

    /** \brief return enriched displacements at the nodes of the element*/
    const Vector & getEnrichedDisplacements() const;

    /** \brief Return the set of eigenvector forming the reference frame of the principal stresses*/
    std::vector<Point> getPrincipalFrame ( const Point &p, bool local = false, VirtualMachine * vm = nullptr, StressCalculationMethod m = REAL_STRESS )  ;

    /** \brief return the linear interpolating factors for the displacement field at the given point*/
    std::vector<double> getInterpolatingFactors ( const Point & p, bool local = false ) const ;

    /** \brief return the linear enrichment interpolating factors for the displacement field at the given point*/
    std::vector<double> getEnrichedInterpolatingFactors ( const Point & p, bool local = false ) const ;
    
    virtual void step ( double dt, const Vector* d ) ;
    virtual bool prepare(const Vector &extrapolatedDisplacements) ;

    double getTime() const ;
    double getDeltaTime() const ;
    double getNodalCentralTime() const ;
    double getNodalDeltaTime() const ;

    IntegrableEntity * getParent() const {
        return parent ;
    }

    virtual void initialize ( Mesh<DelaunayTetrahedron,DelaunayTreeItem3D> * msh ) ;
    virtual void initialize ( Mesh< DelaunayTriangle, DelaunayTreeItem >* msh ) ;

    const Mesh< DelaunayTriangle, DelaunayTreeItem > * getMesh2D() const {
        return mesh2d ;
    } ;
    const Mesh< DelaunayTetrahedron, DelaunayTreeItem3D > * getMesh3D() const {
        return  mesh3d ;
    } ;

    Mesh< DelaunayTriangle, DelaunayTreeItem > * getMesh2D() {
        return mesh2d ;
    } ;
    Mesh< DelaunayTetrahedron, DelaunayTreeItem3D > * getMesh3D() {
        return  mesh3d ;
    } ;

} ;

class ElementStateWithInternalVariables : public ElementState
{
protected:
    std::vector<std::vector<Vector> > internalVariablesAtGaussPoints ;
    int n ;
    int p ;

public:
    ElementStateWithInternalVariables ( IntegrableEntity *, int n, int p ) ;

    int numberOfInternalVariables() const {
        return n ;
    }
    void setNumberOfGaussPoints ( int g ) ;

    int sizeOfInternalVariable() const {
        return p ;
    }

    virtual void getField ( FieldType f, const Point & p, Vector & ret, bool local, VirtualMachine * vm = nullptr, int i = 0 )  ;

    virtual void getField ( FieldType f, const PointArray & p, Vector & ret, bool local, VirtualMachine * vm = nullptr, int i = 0 )   ;

    virtual void getField ( FieldType f, const std::valarray<std::pair<Point, double> > & p, Vector & ret, bool local, VirtualMachine * vm = nullptr, int i = 0 )   ;

// 	virtual void getFieldAtNodes( FieldType f, Vector & ret, int i = 0) ;

    virtual void getFieldAtGaussPoint ( FieldType f, size_t g, Vector & ret, VirtualMachine * vm = nullptr, int i = 0 ) ;

    virtual void getField ( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm = nullptr, int i = 0, int j = 0 )   ;

    virtual void getField ( FieldType f1, FieldType f2, const PointArray & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm = nullptr, int i = 0, int j = 0 )   ;

    virtual void getField ( FieldType f1, FieldType f2, const std::valarray<std::pair<Point, double> > & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm = nullptr, int i = 0, int j = 0 )   ;

// 	virtual void getFieldAtNodes( FieldType f1, FieldType f2, Vector & ret1, Vector & ret2, int i = 0, int j = 0) ;

    virtual void getFieldAtGaussPoint ( FieldType f1, FieldType f2, size_t g, Vector & ret1, Vector & ret2, VirtualMachine * vm = nullptr, int i = 0, int j = 0 ) ;

    virtual void initialize ( Mesh<DelaunayTetrahedron,DelaunayTreeItem3D> * msh ) ;
    virtual void initialize ( Mesh< DelaunayTriangle, DelaunayTreeItem >* msh ) ;

    virtual void setInternalVariableAtGaussPoint ( Vector & v, size_t g, int i ) ;

} ;

/** \brief container for a set of Gauss points*/
struct GaussPointArray {
    std::valarray< std::pair<Point, double> > gaussPoints ;

    bool regularGrid = false ;

    GaussPointArray() : gaussPoints ( std::make_pair ( Point(), 1. ),1 ) { } ;
    GaussPointArray ( const std::pair<Point, double> & p ) : gaussPoints ( p, 1 ) { } ;
    GaussPointArray ( const GaussPointArray & gp ) : gaussPoints ( gp.gaussPoints ) { } ;
    GaussPointArray ( const std::valarray< std::pair<Point, double> > & array/*, int i*/ ) : gaussPoints ( array ) { } ;
    GaussPointArray ( const std::vector< std::pair<Point, double> > & array ) : gaussPoints ( array.size() ) { std::copy(array.begin(), array.end(), &gaussPoints[0]) ; } ;
    void operator = ( const GaussPointArray & gp ) {
        gaussPoints.resize ( gp.gaussPoints.size() ) ;
        gaussPoints = gp.gaussPoints ;
//         id = gp.getId() ;
    }
    
   void operator = ( const std::valarray< std::pair<Point, double> > & gp ) {
        gaussPoints.resize ( gp.size() ) ;
        gaussPoints = gp ;
    }
    
    void operator = ( const std::vector< std::pair<Point, double> > & gp ) {
        gaussPoints.resize ( gp.size() ) ;
        std::copy(gp.begin(), gp.end(), &gaussPoints[0]) ;
    }
} ;

/** \brief Abstract class for the representation of elements
 */
struct IntegrableEntity : public Geometry
{

    bool largeDeformations = false ;
    Order order ;
    ElementState * state ;
    std::vector<BoundaryCondition *> * boundaryConditionCache ;
    GaussPointArray * cachedGps ;
    const GaussPointArray * getCachedGaussPoints() const {
        return cachedGps ;
    };
    void setCachedGaussPoints ( GaussPointArray * gp ) ;

    bool enrichmentUpdated ;
    bool behaviourUpdated ;
    bool behaviourForcesUpdated ;
    bool behaviourViscoUpdated ;
    bool needAssembly ;

    IntegrableEntity() ;
    virtual void getSecondJacobianMatrix ( const Point &p, Matrix & t1, Matrix & t2 ) {
        std::cout << "not implemented..." << std::endl ;
        exit ( 0 ) ;
    } ;
    virtual void getThirdJacobianMatrix ( const Point &p, Matrix & t1, Matrix & t2, Matrix & t3 ) {
        std::cout << "not implemented..." << std::endl ;
        exit ( 0 ) ;
    } ;
    virtual ~IntegrableEntity() ;
    virtual const Point & getPoint ( size_t i ) const = 0 ;
    virtual Point & getPoint ( size_t i )  = 0 ;
    virtual const GaussPointArray & getGaussPoints() = 0;
    virtual Point inLocalCoordinates ( const Point & p ) const  = 0;
    virtual double area() const {
        return 0 ;
    }
    virtual double volume() const {
        return 0 ;
    }

    virtual void addBoundaryCondition ( BoundaryCondition * bc ) ;
    virtual void clearBoundaryConditions() ;
    virtual bool hasBoundaryConditions() {
        return boundaryConditionCache ;
    }

    virtual Function getXTransform() const = 0;
    virtual Function getYTransform() const = 0;
    virtual Function getZTransform() const ;
    virtual Function getTTransform() const ;
    virtual Function getXTransformAtCentralNodalTime() const ;// { return getXTransform() ; }
    virtual Function getYTransformAtCentralNodalTime() const ;//{ return getYTransform() ; }
    virtual Function getZTransformAtCentralNodalTime() const ;//{ return getZTransform() ; }
    virtual Function getTTransformAtCentralNodalTime() const ;//{ return getTTransform() ; }

    virtual const std::valarray< Function >  & getShapeFunctions() const = 0 ;
    virtual std::valarray< Function >  & getShapeFunctions() = 0 ;
    virtual Function & getShapeFunction ( size_t i )  = 0 ;
    virtual const Function & getShapeFunction ( size_t i ) const = 0 ;

    virtual const std::vector< Function> & getEnrichmentFunctions() const = 0 ;
    virtual const Function & getEnrichmentFunction ( size_t i ) const = 0;

    virtual Order getOrder() const  = 0 ;

    virtual std::vector<size_t> clearEnrichment ( const Geometry * g ) = 0 ;
    virtual std::vector<size_t> clearAllEnrichment() = 0;
    virtual const std::vector< size_t > getDofIds() const = 0;

    virtual Form * getBehaviour() const = 0;    
    virtual bool matrixUpdated() const = 0 ;
    
    virtual std::valarray<std::valarray<Matrix> > & getElementaryMatrix( VirtualMachine * vm = nullptr)  = 0 ;
    virtual std::valarray<std::valarray<Matrix> > & getViscousElementaryMatrix(VirtualMachine * vm = nullptr)  = 0 ;
// 	virtual Vector getForces() = 0 ;
    virtual void applyBoundaryCondition ( Assembly * a ) ;

    virtual void adjustElementaryMatrix ( double previousTimeStep, double nextTimeStep ) { } ;

    virtual bool isMoved() const = 0 ;
    virtual void print() const = 0;

    virtual const ElementState & getState() const ;
    virtual ElementState & getState() final;
    virtual void setState ( ElementState * state ) ;

    virtual Mesh<DelaunayTriangle, DelaunayTreeItem> * get2DMesh() const = 0 ;
    virtual Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * get3DMesh() const = 0 ;
    virtual void setOrder ( Order o ) = 0;

} ;

/** \brief A Form for DIM degrees of freedom.
 */
class Form
{
    friend BoundaryCondition ;
protected:
    bool time_d ;
    bool space_d ;
    bool symmetric ;
    size_t num_dof ;

    const Geometry * source ;
    CollisionDetector * collisionDetection = nullptr;
    GeometryBasedEffect * contactModel = nullptr;
    

public:

    /** A form has at least a parameter, which takes the shape of a Matrix*/
    Matrix param ;
    /** The type helps to know the available parameters and methods of the subclasses*/
    ParametersType type;

    Form ( const Matrix & p, bool t = false, bool s = false, size_t numdof = 2, bool sym = true ) : time_d ( t ), space_d ( s ), symmetric ( sym ), num_dof ( numdof ),  source ( nullptr ), param ( p ) { } ;
    Form() : time_d ( false ), space_d ( false ), symmetric ( false ), num_dof ( 0 ), source ( nullptr ), param ( Matrix ( 0,0 ) ) { } ;

    /** apply the form on a pair of functions
     *
     * @param p_i First form function
     * @param p_j Second Form Function
     * @param e The element in which the integration is performed
     * @return The vector of the symbolic matrices at the integration points
     */
    virtual void apply ( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &, VirtualMachine * vm ) const = 0 ;
    virtual void applyViscous ( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &, VirtualMachine * vm ) const { } ;

    virtual bool isSymmetric() const final{
        return symmetric ;
    }
    virtual void setSymmetric ( bool s ) final{
        symmetric = s ;
    }

    virtual bool timeDependent() const final{
        return this->time_d ;
    } ;

    virtual void setTimeDependent ( bool t ) final{
        time_d = t ;
    }
    virtual void setSpaceDependent ( bool s ) final{
        space_d = s ;
    }

    virtual bool spaceDependent() const final{
        return this->space_d ;
    } ;

    virtual Vector getForcesFromAppliedStress ( const Vector & data, const Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, const std::vector<Variable> & v, bool isVolumic = false, const Vector & normal = Vector(), VirtualMachine *vm = nullptr ) const ;

    virtual Vector getForcesFromAppliedStress ( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic = false, const Vector & normal = Vector() , VirtualMachine *vm = nullptr) ;

    virtual const Geometry * getSource() const final {
        return source  ;
    }
    virtual void setSource ( const Geometry *  src ) final;

    virtual bool hasInducedForces() const ;

    virtual bool isViscous() const {
        return false ;
    }

    virtual size_t getNumberOfDegreesOfFreedom() const final {
        return num_dof ;
    }

    virtual void transform ( const ElementarySurface * ) { } ;
    virtual void transform ( const ElementaryVolume * ) { } ;

    /** Step through time
     *
     * @param timestep length of the timestep
     * @param currentState State of the element with this behaviour
     */
    virtual void step ( double timestep, ElementState & currentState, double maxScore ) = 0;

    virtual bool fractured() const = 0 ;
    virtual bool changed() const {
        return false ;
    } ;
    virtual void getForces (  ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & v )  { };
    virtual std::vector<BoundaryCondition * > getBoundaryConditions ( const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv ) const ;

    virtual Form * getCopy() const = 0 ;

    virtual Matrix getTensor ( const Point & p, IntegrableEntity * e = nullptr, int g = -1 ) const {
        return param ;
    }

    virtual Matrix getViscousTensor ( const Point & p, IntegrableEntity * e = nullptr, int g = -1 ) const {
        return param * 0 ;
    }

    virtual Matrix getTensorDot ( const Point & p, double dt, bool cacl = false ) const ; 

    virtual Matrix getViscousTensorDot ( const Point & p, double dt, bool cacl = false ) const ; 

    virtual void getTensorDotAtGaussPoints( const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<std::pair<Matrix, Matrix> > & ret, bool calc = false) const ;

    virtual void getViscousTensorDotAtGaussPoints( const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<std::pair<Matrix, Matrix> > & ret, bool calc = false) const ;

    virtual void setTensor ( const Matrix & m ) {
        param = m ;
    }

    virtual Matrix getPreviousTensor ( const Point & p ) const {
        return param ;
    }

    virtual Vector getImposedStress ( const Point & p, IntegrableEntity * e = nullptr, int g = -1 ) const ;
    virtual Vector getImposedStrain ( const Point & p, IntegrableEntity * e = nullptr, int g = -1 ) const ;

    virtual ~Form() {

    } ;

    virtual FractureCriterion * getFractureCriterion() const {
        return nullptr ;
    }
    
    virtual CollisionDetector * getCollisionDetection() const {
        return collisionDetection ;
    }

    virtual void setFractureCriterion ( FractureCriterion * frac ) { }

    virtual DamageModel * getDamageModel() const {
        return nullptr ;
    }
    
    virtual GeometryBasedEffect * getGeometryBasedEffect() const {
        return contactModel ;
    }

    virtual ElementState * createElementState ( IntegrableEntity * e ) ;

    virtual void preProcess ( double timeStep, ElementState & currentState ) ;

} ;

Matrix makeStressOrStrainMatrix ( const Vector & stressOrStrain ) ;

int isGaussPoint ( const Point & p, IntegrableEntity * e ) ;
Vector toPrincipal ( const Vector & stressOrStrain, CompositionType t) ;
} 


#endif
