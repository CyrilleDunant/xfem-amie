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

#ifndef __DIFFUSION_H_
#define __DIFFUSION_H_

#include "physics_base.h"
#include "homogenization/homogenization_base.h"

namespace Amie
{

/** \brief A linear Diffusion Law
* The field param is the diffusion matrix
*/
struct Diffusion : public LinearForm
{
    std::vector<Variable> v ;
    /** \brief Constructor
    *
    * @param rig Complete expression of the diffusion Tensor
    */
    Diffusion(const Matrix & rig) ;

    virtual ~Diffusion() ;

    /** \brief Apply the law.
     *
     * @param p_i first basis function.
     * @param p_j second basis function.
     * @param gp Gauss points to consider
     * @param Jinv Inverse Jacobian Matrices at the Gauss points
     * @param ret store the result
     * @param vm VirtualMachine to use
     */
    virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    virtual void applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    /** \brief return false
     *
     * @return false
     */
    virtual bool fractured() const ;

    virtual bool isViscous() const {
        return true ;
    }

    /** \brief Return a Copy of the behaviour
     *
     * @return a new Diffusion
     */
    virtual Form * getCopy() const ;

} ;

}

#endif
