//
// C++ Interface: stiffness_with_imposed_deformation_and_fracture
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __STIFFNESS_WITH_DEF_FRACT
#define __STIFFNESS_WITH_DEF_FRACT

#include "physics_base.h"
#include "damagemodels/isotropiclineardamage.h"
#include "fracturecriteria/fracturecriterion.h"

namespace Amie
{


/** \brief A linear elastic law with an imposed strain
*
* The field param is the Cauchy-Green Strain Tensor
* The imposed deformation are given as a vector
*/
struct StiffnessWithImposedDeformationAndFracture : public LinearForm
{
    std::vector<Variable> v ;
    Vector imposed ;
    IsotropicLinearDamage * dfunc ;
    double eps ;
    Vector previousPreviousDamage ;
    Vector intermediateDamage ;
    Vector previousDamage ;
    int count ;
    FractureCriterion * criterion ;
    bool frac ;
    bool change ;
    bool previouschange ;
    double sigmaRupt ;
    double init ;
    Vector damage ;

    /** \brief Constructor
    *
    * @param rig Complete expression of the Cauchy-Green Strain Tensor
    * @param imposedDef Imposed deformation
    */
    StiffnessWithImposedDeformationAndFracture(const Matrix & rig, Vector imposedDef, FractureCriterion * crit) ;

    virtual ~StiffnessWithImposedDeformationAndFracture() ;

    /** Apply the law.
    * @param p_i first basis polynomial.
    * @param p_j second basis polynomial.
    * @param gp integration points
    * @param Jinv inverse Jacobian matrix at the integration points
    * @param ret, matrix to store the result
    * @param vm VirtualMachine to use to perform the computation
    * @return matrix resulting of \f$ \nabla H^T K \nabla H \f$.
    */
    virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    /** \brief return false */
    virtual bool fractured() const ;

    /** \brief Return a copy of this Behaviour*/
    virtual Form * getCopy() const ;

    /** \brief Return the Vector of imposed Stress at the considered point. As the imposed stress is uniform, the point is ignored*/
    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

    /** \brief Acessor, return the stifness criterion in use*/
    virtual FractureCriterion * getFractureCriterion() const ;

    virtual DamageModel * getDamageModel() const ;

    /** \brief Check for fracture
    *
    * @param timestep elapsed time
    * @param currentState state of the element
    *
    * if the yield criterion is true, se fractured state to true
    */
    virtual void step(double timestep, ElementState & currentState, double maxscore) ;

    /** \brief return true if the damage state has been modfied*/
    virtual bool changed() const ;

    /** \brief Return the (damaged) Stifness tensor*/
    virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;




} ;


}

#endif
