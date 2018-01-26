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

#ifndef __STIFFNESS_WITH_LARGE_DEF
#define __STIFFNESS_WITH_LARGE_DEF

#include "physics_base.h"

namespace Amie
{


/** \brief A linear elastic law with an imposed strain
*
* The field param is the Cauchy-Green Strain Tensor
* The imposed deformation are given as a vector
*/

/*PARSE StiffnessWithLargeDeformation Form
    @value[young_modulus] // value of the Young modulus
    @value[poisson_ratio] // value of the Poisson ratio
    @value[imposed_deformation] // value of the linear imposed strain
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
    @string<IsotropicMaterialParameters>[material_parameters] YOUNG_POISSON // describes how to build the stiffness matrix
*/
struct StiffnessWithLargeDeformation : public LinearForm
{
    friend class Triangle ;
    std::vector<Variable> v ;
    Vector imposed ;
    bool active = false ;
    bool change = false ;
    double initialVolume ;
    double poisson ; 
    Matrix largeDeformationTransformc ;
    double globalTransformAngle ;

    std::vector<double> initialPosition ;
    /** \brief Constructor
    *
    * @param rig Complete expression of the Cauchy-Green Strain Tensor
    * @param imposedDef Imposed deformation
    */
    StiffnessWithLargeDeformation(const Matrix & rig) ;

    StiffnessWithLargeDeformation(double E, double nu, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, IsotropicMaterialParameters hooke = YOUNG_POISSON) ;

    virtual ~StiffnessWithLargeDeformation() ;

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
    
    virtual bool changed() const ;

    /** \brief Return a copy of this Behaviour*/
    virtual Form * getCopy() const ;

    virtual bool hasInducedForces() const {
        return true ;
    }

    /** \brief Return the Vector of imposed Stress at the considered point. As the imposed stress is uniform, the point is ignored*/
    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e, int g = -1) const ;
    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e, int g = -1) const ;

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

    virtual void step(double timestep, ElementState & currentState, double maxScore) ;

} ;


}

#endif
