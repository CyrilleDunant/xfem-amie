//
// C++ Interface: isotropiclineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MU_PLASTIC_STRAIN_H
#define MU_PLASTIC_STRAIN_H

#include "damagemodel.h"

namespace Mu {

/** \brief Isotropic linear damage model. The stifness of an affected element is scaled by a factor between 1 and 1 - .9999999
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class PlasticStrain : public DamageModel
{
protected:
	ElementState * es ;
	std::vector<Variable> v ;
	Matrix * param ;
	Vector imposedStrain ;
	Vector previousImposedStrain ;
	Vector lastStress ;
	double c_psi ;
	double plasticVariable ; 
	double eps_f ;
	double plasticFlowPotential(const Matrix & m) const ;
	
public:
	/** \brief Constructor. Set the number of degrees of freedom
	 * 
	 * @param numDof number of degrees of freedom
	 */
	PlasticStrain() ;

	virtual ~PlasticStrain();

	/** \brief Increment the damage
	 * 
	 * @param s ElementState passed as a parameter
	 */
	virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;

	/** \brief compute the new stifness matrix after damage
	 * 
	 * @param m Matrix to modify
	 * @return m
	 */
	virtual Matrix apply(const Matrix & m) const;
	
	virtual bool hasInducedBoundaryConditions() const {return true ;} ;
	virtual bool hasInducedForces() const {return true ;}
	
	virtual Matrix applyPrevious(const Matrix & m) const;

	/** \brief return true is the element concerned is fractured 
		*/
	virtual bool fractured() const  ;
	
// 	virtual void step(ElementState & s) ;
	
	virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
	virtual Vector getImposedStress(const Point & p) const ;
	
	virtual DamageModel * getCopy() const { return new PlasticStrain() ;}
	
	virtual void postProcess() ;
	double getDamage() const ;
};

}

#endif
