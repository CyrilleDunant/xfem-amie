//
// C++ Interface: stiffness_with_imposed_deformation
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __STIFFNESS_WITH_VARIABLE_DEF
#define __STIFFNESS_WITH_VARIABLE_DEF

#include "physics_base.h"

namespace Amie
{


/** \brief A linear Elastic Law
* The field param is the Cauchy-Green Strain Tensor
* The imposed deformation are given as a vector
* They are incremented as a fraction of the timestep
*/
struct StiffnessWithVariableImposedDeformation : public LinearForm
{
    Vector imposed ;
    Function x ;
    Function y ;
    std::vector<Variable> v ;

    /** Constructor
    *
    * @param rig Complete expression of the Cauchy-Green Strain Tensor
    */
    StiffnessWithVariableImposedDeformation(const Matrix & rig, Vector imposedDef) ;

    virtual ~StiffnessWithVariableImposedDeformation() ;

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

    /** \brief Return false*/
    virtual bool fractured() const ;

    /** \brief Return a copy of this Behaviour. The current imposed strain is also copied*/
    virtual Form * getCopy() const ;

    /** \brief Return the current imposed strain vector*/
    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    /** \brief Increment the loading. The loading is increased by the timestep*a random value between 0 and 1 .*/
    virtual void step(double timestep, ElementState & currentState, double maxscore) ;


    std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

} ;


}

#endif
