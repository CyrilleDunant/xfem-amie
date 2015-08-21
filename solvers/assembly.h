
//
// C++ Interface: assembly
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __ASSEMBLY_NEW_H_
#define __ASSEMBLY_NEW_H_

#include <map>
#include <set>
#include <complex>
#include <iomanip>

#include "../mesher/delaunay.h"
#include "../mesher/delaunay_3d.h"
#include "../elements/elements.h"
#include "../physics/physics_base.h"
#include "../sparse/sparse_matrix.h"


namespace Amie
{

class LinearSolver ;
class Preconditionner ;


typedef enum
{
    GENERAL = 0,
    FIX_ALONG_ALL,
    FIX_ALONG_XI,
    SET_ALONG_XI,
    INCREMENT_ALONG_XI,
    FIX_ALONG_ETA,
    SET_ALONG_ETA,
    INCREMENT_ALONG_ETA,
    FIX_ALONG_ZETA,
    SET_ALONG_ZETA,
    INCREMENT_ALONG_ZETA,
    FIX_ALONG_XI_ETA,
    SET_ALONG_XI_ETA,
    INCREMENT_ALONG_XI_ETA,
    FIX_ALONG_XI_ZETA,
    SET_ALONG_XI_ZETA,
    INCREMENT_ALONG_XI_ZETA,
    FIX_ALONG_ETA_ZETA,
    SET_ALONG_ETA_ZETA,
    INCREMENT_ALONG_ETA_ZETA,
    FIX_ALONG_INDEXED_AXIS,
    SET_ALONG_INDEXED_AXIS,
    INCREMENT_ALONG_INDEXED_AXIS,
    SET_PROPORTIONAL_DISPLACEMENT,
    SET_PROPORTIONAL_DISPLACEMENT_XI_ETA,
    SET_PROPORTIONAL_DISPLACEMENT_XI_ZETA,
    SET_PROPORTIONAL_DISPLACEMENT_ETA_XI,
    SET_PROPORTIONAL_DISPLACEMENT_ETA_ZETA,
    SET_PROPORTIONAL_DISPLACEMENT_ZETA_XI,
    SET_PROPORTIONAL_DISPLACEMENT_ZETA_ETA,
    FIX_NORMAL_DISPLACEMENT,
    FIX_TANGENT_DISPLACEMENT,
    SET_NORMAL_DISPLACEMENT,
    SET_TANGENT_DISPLACEMENT,
    SET_FORCE_XI,
    SET_FORCE_ETA,
    SET_FORCE_ZETA,
    SET_FORCE_INDEXED_AXIS,
    SET_FLUX_XI,
    SET_FLUX_ETA,
    SET_FLUX_ZETA,
    SET_VOLUMIC_STRESS_XI,
    SET_VOLUMIC_STRESS_XI_ETA,
    SET_VOLUMIC_STRESS_ETA,
    SET_VOLUMIC_STRESS_ZETA,
    SET_STRESS_XI,
    SET_STRESS_ETA,
    SET_STRESS_ZETA,
    SET_NORMAL_STRESS,
    SET_TANGENT_STRESS,
    VERTICAL_PLANE_SECTIONS,
    HORIZONTAL_PLANE_SECTIONS,
    nullptr_CONDITION,
    SET_GLOBAL_FORCE_VECTOR,
} LagrangeMultiplierType ;

/** \brief Abstract representation of a Boundary condition. Can be an actual Lagrange Multiplier, or a hint for set displacements or forces.*/
class LagrangeMultiplier
{
public:
    std::valarray<unsigned int> ids ;
    Vector coefs ;
    int id ;
    double value ;
    std::vector<std::pair<unsigned int, double> > hints ;
public:

    /** \brief Constructor. Create a Lagrange Multiplier for a set of dofs, related by a linear form.
    *
    * @param i Indices of the dofs
    * @param c coefficients associated to the dofs
    * @param v right-hand side of the equation
    * @param my_id this LM's ID
    */
    LagrangeMultiplier(std::valarray<unsigned int> i, Vector c, double v ,int my_id  = -1) ;

    /** \brief Create a dummy LM with given ID*/
    LagrangeMultiplier(int _id) {
        id  = _id ;
    }

    /** \brief Return the sparse matrix corresponding to this lagrange Multiplier*/
    CoordinateIndexedIncompleteSparseMatrix getMatrix() const ;

    /** \brief set the id*/
    void setId(const int) ;

    /** \brief get the ID*/
    int getId() const ;

    /** \brief get all affected dofs*/
    const std::valarray<unsigned int> getDofIds() const ;

    /** \brief get Right-hand-side value*/
    double getValue() const  ;

    /** \brief get linear form*/
    std::vector<std::pair<unsigned int, double> > getHints() const ;

    /** \brief set linear form*/
    void setHint(std::pair<unsigned int, double>) ;

    bool operator < (const LagrangeMultiplier & m) const
    {
        return id < m.getId() ;
    }

    bool operator < (const int & m) const
    {
        return id < m ;
    }
    
    bool operator > (const LagrangeMultiplier & m) const
    {
        return id > m.getId() ;
    }

    bool operator > (const int & m) const
    {
        return id > m ;
    }

    bool operator == (const LagrangeMultiplier & m) const
    {
        return id == m.getId() ;
    }

    bool operator == (const int & m) const
    {
        return id == m ;
    }

    void print() const
    {
        std::cout << "Multiplier affectig ID " << id << ", val = " << value << std::endl ;
    }

    LagrangeMultiplierType type ;
} ;

bool operator < (const int i, const LagrangeMultiplier & l) ;

/** \brief functor for usage with STL containers. order LagrangeMultipliers*/
struct MultiplierHasIdSupEq
{
    int m_id ;
    MultiplierHasIdSupEq(int id) {
        m_id = id ;
    }
    bool operator()(const LagrangeMultiplier & m)
    {
        return m_id >= m.getId() ;
    }
} ;

/** \brief functor for usage with STL containers. order LagrangeMultipliers*/
struct MultiplierHasIdInfEq
{
    int m_id ;
    MultiplierHasIdInfEq(int id) {
        m_id = id ;
    }
    bool operator()(const LagrangeMultiplier & m)
    {
        return m_id <= m.getId() ;
    }
} ;

/** \brief functor for usage with STL containers. order LagrangeMultipliers*/
struct MultiplierHasId
{
    int m_id ;
    MultiplierHasId(int id) {
        m_id = id ;
    }
    bool operator()(const LagrangeMultiplier & m)
    {
        return m_id == m.getId() ;
    }
} ;

/** \brief Assembly of the elementary Matrices
*
* Elements are add() ed to the assembly. The assembly provides methods for the setting of boundary conditions.
* A solver is finally called.
*/
class Assembly
{
protected:

    std::vector<ElementarySurface *> element2d ;
    std::vector<ElementaryVolume *> element3d ;
    std::vector<double> scales ;
    std::vector<LagrangeMultiplier> multipliers ;
    std::vector<LagrangeMultiplier> multipliersBuffer ;
    size_t ndof ;
    size_t ndofmax ;
    size_t colstart ;
    size_t rowstart ;
    SpaceDimensionality dim ;
    bool removeZeroOnlyLines = true ;

    Vector naturalBoundaryConditionForces ;
    Vector displacements ;
    Vector prevDisplacements ;
    Vector realPrevDisplacements ;
    Vector externalForces ;
    Vector nonLinearExternalForces ;
    Vector addToExternalForces ;
    CoordinateIndexedSparseMatrix * coordinateIndexedMatrix ;
    CoordinateIndexedSparseMaskMatrix * mask ;
    CoordinateIndexedIncompleteSparseMatrix * nonLinearPartialMatrix ;
    std::map<std::pair<size_t, size_t>, double > * boundaryMatrix ;

    bool make_final() ;

    size_t multiplier_offset ;

    double epsilon ;

public:
    Assembly() ;

    virtual ~Assembly() ;

    size_t getMaxDofID() const ;
    size_t getMaxNodeID() const ;
    std::set<unsigned long int> updatedDofs ;

    /** \brief add element to assembly*/
    void add(ElementarySurface * e, double scale = 1.) ;

    /** \brief add element to assembly*/
    void add(ElementaryVolume * e, double scale = 1.) ;

    /** \brief chek if a line is equal to 0 or not*/
    void checkZeroLines() ;

    void setMaxDof(size_t n) {
        ndofmax = n ;
    }

    void setRemoveZeroOnlyLines(bool r) {
        removeZeroOnlyLines = r ;
    }

    /** \brief return true if Assembly is made of 2D elements*/
    bool has2DElements() const {
        return element2d.size() > 0 ;
    }

    /** \brief return true if Assembly is made of 3D elements*/
    bool has3DElements() const {
        return element3d.size() > 0 ;
    }

    /** \brief print assembled matrix ans vector*/
    void print() ;

    /** \brief print the diagonal of the assembled matrix*/
    void printDiag() const ;

    /** \brief apply Boundary conditions*/
    void setBoundaryConditions() ;

    /** \brief perform a sub-step in a non-linear problem*/
    bool nonLinearStep() ;
// 	Vector solve(size_t maxit = 1000) ;

    /** \brief Solve linear system*/
    bool solve(Vector x, size_t maxit = 1000, const bool verbose = false) ;

    /** \brief Solve linear system using Preconditionned Conjugate Gradient (linear/non linear/biconjugate is automatically selected)*/
    bool cgsolve(Vector x0 = Vector(0), int maxit = -1, bool verbose = true) ;

    void setEpsilon(double e) {
        epsilon = e ;
    }
    double getEpsilon() const {
        return epsilon ;
    }

    /** \brief Solve linear system using provided solver*/

    bool mgprepare() ;
    bool mgsolve(LinearSolver * mg, Vector x0 = Vector(0), Preconditionner * pg = nullptr, int maxit = -1) ;

    /** \brief Solve linear system using Conjugate Gradient (linear/non linear/biconjugate is automatically selected)*/
    bool cgnpsolve(Vector b, size_t maxit) ;

    /** \brief return assembled matrix*/
    CoordinateIndexedSparseMatrix & getMatrix() ;

    /** \brief return assembled matrix*/
    const CoordinateIndexedSparseMatrix & getMatrix() const;

    /** \brief return non-linear part of the assembled matrix*/
    CoordinateIndexedIncompleteSparseMatrix & getNonLinearMatrix() ;

    /** \brief return right-and-side vector*/
    Vector & getForces() ;

    /** \brief return virtual forces which come from imposed displacements BCs*/
    Vector & getNaturalBoundaryConditionForces() ;

    /** \brief return right-and-side vector*/
    Vector & getNonLinearForces() ;

    /** \brief return result of the system*/
    Vector & getDisplacements() {
        return displacements ;
    }    
// 	const Vector & getForces() const;

    Vector & getPreviousDisplacements() {
        return realPrevDisplacements ;
    }

    void setPreviousDisplacements() {
        realPrevDisplacements = displacements ;
    }

    /** \brief set boundary condition. Point with ID id is fixed*/
    void fixPoint(size_t id) ;

    /** \brief set boundary condition. Point with ID id has fixed displacement.
    	 * @param id
    	 */
    void setPoint(double ex, size_t id) ;

    /** \brief set boundary condition. Point with ID id has fixed displacement.
    	 * @param ex
    	 * @param ey
    	 * @param id
    	 */
    void setPoint(double ex, double ey, size_t id) ;

    /** \brief set boundary condition. Point with ID id has fixed displacement.
    	 * @param ex
    	 * @param ey
    	 * @param ez
    	 * @param id
    	 */
    void setPoint(double ex, double ey, double ez, size_t id) ;

    /** Set unknown of a given node in a given direction to a given value
    	 * @param var Variable (XI, ETA, or ZETA) along which displacement is fixed
    	 * @param id ID of the node
    	 * @param val value at which unknown is set
     */
    void setPointAlong(Variable, double val, size_t id) ;

    void setPointAlongIndexedAxis(int index, double val, size_t id, bool force = false) ;

    void setPointProportional(Variable, Variable, double val, double offset, size_t id) ;

    void setPointProportional(Variable, std::vector< std::pair< Variable, double > > val, double offset, size_t id) ;

    /** Apply a force on given node
     * @param var Variable (XI, ETA, or ZETA) along which force is applied
     * @param id ID of the node
     * @param val intensity of the force
     */
    void setForceOn(Variable var, double val, size_t id) ;

    /** Apply a force on given node. The force specified is added to the already applied force (if any).
     * @param var Variable (XI, ETA, or ZETA) along which force is applied
     * @param id ID of the node
     * @param val intensity of the force
     */
    void addForceOn(Variable var, double val, size_t id) ;

    void addForceOnIndexedAxis(int index, double val, size_t id) ;

    void addForceToExternalForces( int axis, double val, size_t id ) ;

    void setDisplacementByDof(size_t dof, double val) ;

    /** Apply a force on given node (1D)
     * @param id ID of the node
     * @param val intensity of the force
     */
    void setForceOn(double val, size_t id) ;

    /**Add a Lagrange Multiplier to the problem
     * @param l Multiplier to add
     */
    void addMultiplier(const LagrangeMultiplier & l) ;

    void addForceVector(const Vector & v) ;

    /** \brief set boundary condition. Point with ID id has 0 displacement along prescribed axis*/
    void fixPoint(size_t id, Amie::Variable v) ;

    /** \brief The two points are the same*/
    void setPeriodicPoint(size_t id0, size_t id1) ;

    /** \brief return Froebenius norm of the assembled matrix*/
    double froebeniusNorm() ;

    void setNumberOfDregreesOfFreedom(int dof) ;
    int getNumberOfDegreesOfFreedom() const ;

    void setSpaceDimension(SpaceDimensionality d) ;
    SpaceDimensionality getSpaceDimension() const ;

    void clear() ;
    void clearElements() ;
} ;

/** \brief Assembly of the elementary Matrices
*
* Elements are add() ed to the assembly. The assembly provides methods for the setting of boundary conditions.
* A solver is finally called.
*/
class ParallelAssembly
{
protected:

    std::vector<Assembly> assembly ;
    std::vector<Geometry *> domains ;

public:
    ParallelAssembly(const std::vector<Geometry *> & domains) ;

    virtual ~ParallelAssembly() ;

    size_t getMaxDofID() const ;

    /** \brief add element to assembly*/
    void add(ElementarySurface * e, double scale = 1.) ;

    /** \brief add element to assembly*/
    void add(ElementaryVolume * e, double scale = 1.) ;

    void setMaxDof(size_t n) ;

    void setRemoveZeroOnlyLines(bool r) ;

    /** \brief return true if Assembly is made of 2D elements*/
    bool has2DElements() const ;

    /** \brief return true if Assembly is made of 3D elements*/
    bool has3DElements() const ;

    /** \brief print assembled matrix ans vector*/
    void print() ;

    /** \brief print the diagonal of the assembled matrix*/
    void printDiag() const ;

    /** \brief apply Boundary conditions*/
    void setBoundaryConditions() ;

    /** \brief perform a sub-step in a non-linear problem*/
    bool nonLinearStep() ;
//  Vector solve(size_t maxit = 1000) ;

    /** \brief Solve linear system*/
    bool solve(Vector x, size_t maxit = 1000, const bool verbose = false) ;

    /** \brief Solve linear system using Preconditionned Conjugate Gradient (linear/non linear/biconjugate is automatically selected)*/
    bool cgsolve(Vector x0 = Vector(0), int maxit = -1, bool verbose = true) ;

    void setEpsilon(double e) ;
    double getEpsilon() const ;

    /** \brief Solve linear system using provided solver*/

    bool mgprepare() ;
    bool mgsolve(LinearSolver * mg, Vector x0 = Vector(0), Preconditionner * pg = nullptr, int maxit = -1) ;

    /** \brief Solve linear system using Conjugate Gradient (linear/non linear/biconjugate is automatically selected)*/
    bool cgnpsolve(Vector b, size_t maxit) ;

    /** \brief return assembled matrix*/
    CoordinateIndexedSparseMatrix & getMatrix(int i) ;

    /** \brief return assembled matrix*/
    const CoordinateIndexedSparseMatrix & getMatrix(int i) const;

    /** \brief return non-linear part of the assembled matrix*/
    CoordinateIndexedIncompleteSparseMatrix & getNonLinearMatrix(int i) ;

    /** \brief return right-and-side vector*/
    Vector & getForces(int i) ;

    /** \brief return virtual forces which come from imposed displacements BCs*/
    Vector & getNaturalBoundaryConditionForces(int i) ;

    /** \brief return right-and-side vector*/
    Vector & getNonLinearForces(int i) ;

    /** \brief return result of the system*/
    Vector & getDisplacements(int i) ;
//  const Vector & getForces() const;

    /** \brief set boundary condition. Point with ID id is fixed*/
    void fixPoint(size_t id) ;

    /** \brief set boundary condition. Point with ID id has fixed displacement.
         * @param id
         */
    void setPoint(double ex, size_t id) ;

    /** \brief set boundary condition. Point with ID id has fixed displacement.
         * @param ex
         * @param ey
         * @param id
         */
    void setPoint(double ex, double ey, size_t id) ;

    /** \brief set boundary condition. Point with ID id has fixed displacement.
         * @param ex
         * @param ey
         * @param ez
         * @param id
         */
    void setPoint(double ex, double ey, double ez, size_t id) ;

    /** Set unknown of a given node in a given direction to a given value
         * @param var Variable (XI, ETA, or ZETA) along which displacement is fixed
         * @param id ID of the node
         * @param val value at which unknown is set
     */
    void setPointAlong(Variable, double val, size_t id) ;

    void setPointAlongIndexedAxis(int index, double val, size_t id) ;

    /** Apply a force on given node
     * @param var Variable (XI, ETA, or ZETA) along which force is applied
     * @param id ID of the node
     * @param val intensity of the force
     */
    void setForceOn(Variable var, double val, size_t id) ;

    /** Apply a force on given node. The force specified is added to the already applied force (if any).
     * @param var Variable (XI, ETA, or ZETA) along which force is applied
     * @param id ID of the node
     * @param val intensity of the force
     */
    void addForceOn(Variable var, double val, size_t id) ;

    void addForceOnIndexedAxis(int index, double val, size_t id) ;

    void setDisplacementByDof(size_t dof, double val) ;

    /** Apply a force on given node (1D)
     * @param id ID of the node
     * @param val intensity of the force
     */
    void setForceOn(double val, size_t id) ;

    /**Add a Lagrange Multiplier to the problem
     * @param l Multiplier to add
     */
    void addMultiplier(const LagrangeMultiplier & l) ;

    void addForceVector(const Vector & v) ;

    /** \brief set boundary condition. Point with ID id has 0 displacement along prescribed axis*/
    void fixPoint(size_t id, Amie::Variable v) ;

    /** \brief The two points are the same*/
    void setPeriodicPoint(size_t id0, size_t id1) ;

    /** \brief return Froebenius norm of the assembled matrix*/
    double froebeniusNorm() ;

    void setNumberOfDregreesOfFreedom(int dof) ;
    int getNumberOfDegreesOfFreedom() const ;

    void setSpaceDimension(SpaceDimensionality d) ;
    SpaceDimensionality getSpaceDimension() const ;

    void clear() ;
    void clearElements() ;
} ;




std::map<std::pair<size_t, size_t>, Amie::Matrix> incompleteCholeskyDecomposition(std::map<std::pair<size_t, size_t>, Amie::Matrix>  & morseMatrix) ;
Vector operator *(const std::map< std::pair<size_t, size_t>, Amie::Matrix > A , const Vector x) ;
Vector operator *(const std::map< std::pair<size_t, size_t>, double > A , const Vector x) ;
std::map< std::pair<size_t, size_t>, Amie::Matrix > operator *(const std::map< std::pair<size_t, size_t>, Amie::Matrix > A , const std::map< std::pair<size_t, size_t>, Amie::Matrix > B) ;
}
#endif


