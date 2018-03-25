//
// C++ Interface: stiffness_and_fracture
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __STIFFNESS_LARGE_STRAINS_AND_FRACTURE_H
#define __STIFFNESS_LARGE_STRAIN_AND_FRACTURE_H

#include "physics_base.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/lineardamage.h"
#include "damagemodels/anisotropicdamage.h"
#include "damagemodels/isotropiclineardamage.h"
#include "damagemodels/nonlocalisotropiclineardamage.h"

namespace Amie
{

/** \brief A linear Elastic Law with a pluggable fracture criterion
*
* The field param is the Cauchy-Green Strain Tensor
*/

/*PARSE StiffnessWithLargeDeformationAndFracture Form 
    @value[young_modulus] // value of the Young modulus
    @value[poisson_ratio] // value of the Poisson ratio
    @object<FractureCriterion>[fracture_criterion] // criterion to calculate the damage
    @object<DamageModel>[damage_model] nullptr // algorithm to compute the damage
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
    @string<IsotropicMaterialParameters>[material_parameters] YOUNG_POISSON // describes how to build the stiffness matrix
 */
struct StiffnessWithLargeDeformationAndFracture : public LinearForm
{
    DamageModel * dfunc ;
        Vector imposed ;
    bool active = false ;
    bool change = false ;
    double initialVolume ;
    std::vector<double> initialPosition ;
    
    FractureCriterion * criterion ;
    std::vector<Variable> v ;

    /** \brief Constructor
    *
    * @param rig Complete expression of the Cauchy-Green Strain Tensor
    * @param c  FractureCriterion to use. The behaviour is responsible for deleting the criterion upon cleanup.
    */
    StiffnessWithLargeDeformationAndFracture(const Matrix & rig, FractureCriterion * c, DamageModel * d = nullptr)  ;

    StiffnessWithLargeDeformationAndFracture(double E, double nu, FractureCriterion * c, DamageModel * d = nullptr, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, IsotropicMaterialParameters hooke = YOUNG_POISSON)  ;

    virtual ~StiffnessWithLargeDeformationAndFracture();

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

    /** \brief Check for fracture
    *
    * @param timestep elapsed time
    * @param currentState state of the element
    *
    * if the yield criterion is true, se fractured state to true
    */
    virtual void step(double timestep, ElementState & currentState, double maxscore) ;

    /** \brief Check for fracture state
    *
    * @return true if the element is fractured
    */
    virtual bool fractured() const ;


    /** \brief get Copy of the behaviour
    *
    * The copy also copies the damage state.
    * @return pointer to the copy. Caller is responsible for cleaning memory
    */
    virtual Form * getCopy() const ;

    /** \brief return true if the damage state has been modfied*/
    virtual bool changed() const ;

    /** \brief Return the (damaged) Stiffness tensor*/
    virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    /** \brief Acessor, return the stiffness criterion in use*/
    virtual FractureCriterion * getFractureCriterion() const ;

    virtual DamageModel * getDamageModel() const ;
    
        /** \brief Return the Vector of imposed Stress at the considered point. As the imposed stress is uniform, the point is ignored*/
    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e, int g = -1) const ;
    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e, int g = -1) const ;

    virtual void setFractureCriterion(FractureCriterion * frac) ;
    
    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;



} ;



}


#endif
