// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef __DUAL_BEHAVIOUR_H_
#define __DUAL_BEHAVIOUR_H_

#include "physics_base.h"

namespace Amie
{

/** \brief A Geometry determined dual behaviour.
 * Each Gauss point used is attributed a behaviour depending on whether it lies in or out a given Geometry
 */
class BimaterialInterface final: public LinearForm
{
public:
    Geometry * inGeometry ;
    Form * inBehaviour ;
    Form * outBehaviour ;

    Function xtransform ;
    Function ytransform ;
    Function ztransform ;
    Function ttransform ;

    /** \brief Constructor, set the Behaviour s and the delimiting Geometry
     *
     * @param in Geometry
     * @param inbehaviour
     * @param outbehaviour
     */
    BimaterialInterface(Geometry * in, Form * inbehaviour, Form * outbehaviour);

    virtual ~BimaterialInterface();

    /** \brief Set the coordinate transformation functions to use to determin whether a point expressed in the local coordinates of an element lies in or out the Geometry
     *
     * @param x x transformation
     * @param y y transformation
     */
    virtual void transform(const ElementarySurface*) ;

    /** \brief Set the coordinate transformation functions to use to determin whether a point expressed in the local coordinates of an element lies in or out the Geometry
     *
     * @param x x transformation
     * @param y y transformation
     * @param z z transformation
     */
    virtual void transform(const ElementaryVolume *) final;

    /** \brief return the linear factor of the behaviour corresponding to the position given
     *
     * @param p Point at which to compute the linear parameter
     * @return
     */
    virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Matrix getViscousTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    /** \brief Return the imposed stress at the point considered
     *
     * @param p Point
     * @return stress Vector
     */
    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    /** \brief Apply the behaviour
    	* This overloaded apply() is more efficient and is designed to minimise allocating and dealocating memory.
    	* The result of the computation depends on the location of the Gauss points
    	*
    	* @param p_i first shape function.
    	* @param p_j second shape function.
    	* @param gp Set of gauss points for numerical integration
    	* @param Jinv inverse jacobian matrices at the gauss points
    	* @param ret matrix in which to sore the results
    	* @param vm pointer to the virtual machine dedicated for the computation
    	*/
    virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    virtual void applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    /** \brief get a copy of the behaviour
     * This will create a new Bimateral behaviour, with a copy of both the members.
     * @return pointer to the copy. Caller is responsible fior cleaning memory
     */
    virtual Form * getCopy() const ;

    /** \brief Time-step the behaviour
    * This will step both be behaviours
    * @param timestep delta-time of the step.
    * @param currentState State of the element in which the behaviour is time-stepped
    */
    virtual void step(double timestep, ElementState & currentState, double maxScore) ;

    virtual ElementState * createElementState( IntegrableEntity * e) {
        return inBehaviour->createElementState( e ) ;
    }

    virtual bool hasInducedForces() const {
        return inBehaviour->hasInducedForces() || outBehaviour->hasInducedForces() ;
    }

    bool insideGeometry(const Point & p) const {
        return inGeometry->in(p) ;
    }

    virtual bool isViscous() const {
        return inBehaviour->isViscous() || outBehaviour->isViscous() ;
    }

    virtual bool changed() const ;
    
    virtual bool fractured() const ;

    /** \brief Return the vector of induced forces if any of the behaviours induces internal forces. Return an empty vecor otherwise
    *
    * @param s ElementState to consider
    * @param p_i shape function to consider
    * @param gp Gauss points to use for the computation of the quadrature
    * @param Jinv Inverse Jacobian matrix to use at the Gauss Points
    * @param v Vector in which to store the result
    */
    std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const;

    /** \brief Acessor, return the highest fracture citerion*/
    virtual FractureCriterion * getFractureCriterion() const ;
    virtual DamageModel * getDamageModel() const ;

    
    virtual Vector getForcesFromAppliedStress( const Vector & data, const Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, const std::vector<Variable> & v, bool isVolumic = false, const Vector & normal = Vector(), VirtualMachine *vm = nullptr ) const ;

    virtual Vector getForcesFromAppliedStress(const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic = false, const Vector & normal = Vector() , VirtualMachine *vm = nullptr ) ;

} ;

}

#endif
