//
// C++ Interface: stiffness
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __STIFFNESS_H_
#define __STIFFNESS_H_

#include "physics_base.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/damagemodel.h"
#include "fracturecriteria/nonlocalvonmises.h"

namespace Amie
{

/** \brief A linear Elastic Law
* The field param is the Cauchy-Green Strain Tensor
*/
/*PARSE Stiffness Form 
    @value[young_modulus] // value of the Young modulus
    @value[poisson_ratio] // value of the Poisson ratio
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
 */
struct Stiffness : public LinearForm
{
    std::vector<Variable> v ;
    /** \brief Constructor
    *
    * @param rig Complete expression of the Cauchy-Green Strain Tensor
    */
    Stiffness(const Matrix & rig) ;
    Stiffness(double E, double nu, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS) ;

    virtual ~Stiffness() ;

    /** \brief Apply the law.
     *
     * The matrix is computed as: \f$ \nabla^T h_i K \nabla h_j \f$
     * @param p_i first basis polynomial.
     * @param p_j second basis polynomial.
     * @param gp Gauss Points used for the quadrature
     * @param Jinv Inverse Jacobian Matrices corresponding to the gauss points
     * @param ret Matrix to store the result
     * @param vm virtualMachine to use to compute the result
     */
    virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    /** \brief Return false.*/
    virtual bool fractured() const ;

    /** \brief Return a copy of the behaviour*/
    virtual Form * getCopy() const ;

} ;

struct WeibullDistributedElasticStiffness : public Stiffness 
{
    double variability ;
    WeibullDistributedElasticStiffness(const Matrix & rig, double v) : Stiffness(rig), variability(v) { }
    WeibullDistributedElasticStiffness(double E, double nu, double v, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS) : Stiffness(E, nu, dim, pt), variability(v) { } 

    virtual Form * getCopy() const ;

} ;


/** \brief A linear Elastic Law
* The field param is the Cauchy-Green Strain Tensor
*/
struct DerivedStiffness : public LinearForm
{
    const Matrix & rig ;
    std::vector<Variable> v ;
    /** \brief Constructor
    *
    * @param rig Complete expression of the Cauchy-Green Strain Tensor
    */
    DerivedStiffness(const Matrix & rig) ;

    virtual ~DerivedStiffness() ;

    /** \brief Apply the law.
     *
     * The matrix is computed as: \f$ \nabla^T h_i K \nabla h_j \f$
     * @param p_i first basis polynomial.
     * @param p_j second basis polynomial.
     * @param gp Gauss Points used for the quadrature
     * @param Jinv Inverse Jacobian Matrices corresponding to the gauss points
     * @param ret Matrix to store the result
     * @param vm virtualMachine to use to compute the result
     */
    virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    virtual Matrix getTensor ( const Point & p, IntegrableEntity * e = nullptr, int g = -1 ) const {
        return rig ;
    }

    /** \brief Return false.*/
    virtual bool fractured() const ;

    /** \brief Return a copy of the behaviour*/
    virtual Form * getCopy() const ;

} ;


/** \brief A linear pseudo-plastic law.
* The behaviour is updated when step is called and the secant stifness is computed as a function  of the strain
*/
class PseudoPlastic : public LinearForm
{
protected:
    void fixLastDamage() ;
    NonLocalVonMises  * vm ;
    bool initialised ;
    Vector imposedStrain ;
    Vector previousImposedStrain ;
    double stiffness ;
public:
    std::vector<DelaunayTriangle *> cache ;
    std::vector<Variable> v ;
    double alpha ;
    double lastDamage ;
    bool change ;
    bool frac ;
    bool fixedfrac ;
    double radius ;
    double limitStrain ;

    /** \brief Constructor
    *
    * @param rig Complete expression of the Cauchy-Green Strain Tensor
    */
    PseudoPlastic(const Matrix & rig, double E, double limitStrain, double radius) ;

    virtual ~PseudoPlastic() ;

    /** \brief Apply the law.
     *
     * The matrix is computed as: \f$ \nabla^T h_i K \nabla h_j \f$
     * @param p_i first basis polynomial.
     * @param p_j second basis polynomial.
     * @param gp Gauss Points used for the quadrature
     * @param Jinv Inverse Jacobian Matrices corresponding to the gauss points
     * @param ret Matrix to store the result
     * @param vm virtualMachine to use to compute the result
     */
    virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    /** \brief update behaviour
    *
    * @param timestep elapsed time
    * @param currentState state of the element
    *
    */
    virtual void step(double timestep, ElementState & currentState, double maxScore) ;

    virtual bool changed() const ;
    virtual Matrix getTensor(const Point & p, IntegrableEntity * e = 0, int g = -1) const ;
    virtual Matrix getPreviousTensor(const Point & p) const ;
    virtual FractureCriterion * getFractureCriterion() const ;

    /** \brief Return false.*/
    virtual bool fractured() const ;

    /** \brief Return a copy of the behaviour*/
    virtual Form * getCopy() const ;


} ;

}

#endif
