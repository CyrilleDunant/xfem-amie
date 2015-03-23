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

#ifndef __ORTHO_STIFFNESS_H_
#define __ORTHO_STIFFNESS_H_

#include "physics_base.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/damagemodel.h"
#include "homogenization/homogenization_base.h"
#include "fracturecriteria/nonlocalvonmises.h"

namespace Amie
{

/** \brief A linear Elastic Law
* The field param is the Cauchy-Green Strain Tensor
*/
struct OrthotropicStiffness : public LinearForm
{
public:
    Matrix transform ;
    Matrix transformt ;
    Matrix paramBase ;
    bool transformset ;
    std::vector<Variable> v ;
public:
    /** \brief Constructor
     *
     * @param E_1 stifness in the first principal direction
     * @param E_2 stifness in the second principal direction
     * @param nu_12 Poisson ratio
     * @param nu_21 Poisson ratio
     * @param angle angle of the fibres
     * @param poissondefined dummy parameter
     */
    OrthotropicStiffness(double E_1, double E_2, double nu_12,  double nu_21, double angle, bool poissondefined) ;

    OrthotropicStiffness(const OrthotropicStiffness * source) ;

    virtual Matrix getTensor(const Point& p, IntegrableEntity* e = nullptr, int g = -1) const;

    void setAngle(double) ;
    void setStiffness(double E_1, double E_2, double G, double nu) ;
    void setStiffness(double E_1, double E_2, double E_3, double G_1, double G_2, double G_3, double nu) ;
    /** \brief Constructor
    *
    * @param E_1 stifness in the first principal direction
    * @param E_2 stifness in the second principal direction
    * @param G shear modulus
    * @param nu 12 Poisson ratio
    * @param angle angle of the fibres
    */
    OrthotropicStiffness(double E_1, double E_2, double G,  double nu, double angle) ;
    /** \brief Constructor
    *
    * @param E_1 stifness in the first principal direction
    * @param E_2 stifness in the second principal direction
    * @param E_3 stifness in the third principal direction
    * @param G_1 shear modulus
    * @param G_2 shear modulus
    * @param G_3 shear modulus
    * @param nu 12 Poisson ratio
    * @param angle angle of the fibres
    */
    OrthotropicStiffness(double E_1, double E_2, double E_3, double G_1, double G_2, double G_3,  double nu, double angle) ;

    virtual ~OrthotropicStiffness() ;


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
}

#endif
