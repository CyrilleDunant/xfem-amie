//
// C++ Interface: weibull_distributed_stiffness
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __WEIBULL_STIFFNESS_H_
#define __WEIBULL_STIFFNESS_H_

#include "physics_base.h"
#include "damagemodels/damagemodel.h"
#include "fracturecriteria/fracturecriterion.h"

namespace Amie
{

/** \brief A linear Elastic Law
* The field param is the Cauchy-Green Strain Tensor
* Actual tensors are distributed around the prescribed value following a Weibull distribution
*/
/*PARSE WeibullDistributedStiffness Form 
    @value[young_modulus] // value of the Young modulus
    @value[poisson_ratio] // value of the Poisson ratio
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @value[compressive_strength] // value of the compressive strength of the material
    @value[tensile_strength] // value of the tensile strength of the material
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
    @value[variability] 0.2 // variability of the mechanical properties
    @value[material_characteristic_radius] 0.001 // radius of the non-local damage model
    @string<IsotropicMaterialParameters>[material_parameters] YOUNG_POISSON // describes how to build the stiffness matrix
 */
struct WeibullDistributedStiffness : public LinearForm
{
    double materialRadius ;
    std::vector<Variable> v ;
    double variability ;
    double down ;
    double up ;
    double E ;
    double nu ;
    SpaceDimensionality dim ;

    DamageModel * damageModel ;
    /** Constructor
    *
    * @param rig Complete expression of the Cauchy-Green Strain Tensor
    * @param cri stress limit for the Mohr - Coulomb criterion to use
    */
    WeibullDistributedStiffness(double E, double nu, SpaceDimensionality dim, double down, double up, planeType pt = PLANE_STRAIN, double var = 0.2, double radius = 0.001, IsotropicMaterialParameters hooke = YOUNG_POISSON)  ;

    void setNeighbourhoodRadius(double r) {
        materialRadius = r ;
    }

    virtual ~WeibullDistributedStiffness();

    /** \brief Apply the behaviour
    	* This overloaded apply() is more efficient and is designed to minimise allocating and dealocating memory.
    	* This method  is never used, however, because the behaviour only serves as a template for copying
    	*
    	* @param p_i first shape function.
    	* @param p_j second shape function.
    	* @param gp Set of gauss points for numerical integration
    	* @param Jinv inverse jacobian matrices at the gauss points
    	* @param ret matrix in which to sore the results
    	* @param vm pointer to the virtual machine dedicated for the computation
    	*/
    virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    virtual bool fractured() const ;

    /** \brief return a StifnessAndFracture Behaviour
    *
    * The parameters for the copy are determined randomly using a Weibull distribution. As this behaviour is spaceDependant, only the copies are used in elements
    */
    virtual Form * getCopy() const ;

    void setDamageModel(DamageModel * newmodel) ;


} ;

}

#endif
