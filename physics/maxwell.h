//
// C++ Interface: generalized maxwell model with finite difference time-stepping
//
// Description: 
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2012
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __MAXWELL_H_
#define __MAXWELL_H_

#include "physics_base.h"
#include "stiffness.h"



namespace Mu
{
  
struct Maxwell : public LinearForm
{
	Matrix decay ;
	std::vector<Variable> v ;
	/** \brief Constructor
	*
	* @param rig Complete expression of the Cauchy-Green Strain Tensor
	* @param eta Complete expression of the viscosity Tensor
	*/
	Maxwell( const Matrix & rig, const Matrix & eta ) ;

	virtual ~Maxwell() ;

	/** \brief Apply the law.
	 *
	 * The matrix is computed as: \f$ \nabla^T h_i K \nabla h_j + \dot{\nabla}^T h_i E \dot{\nabla} h_j\f$
	 * @param p_i first basis polynomial.
	 * @param p_j second basis polynomial.
	 * @param gp Gauss Points used for the quadrature
	 * @param Jinv Inverse Jacobian Matrices corresponding to the gauss points
	 * @param ret Matrix to store the result
	 * @param vm virtualMachine to use to compute the result
	 */
	virtual void apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const ;

	/** \brief is Element fractured
	 *
	 * @return false
	 */
	virtual bool fractured() const ;

	virtual ElementState * createElementState( IntegrableEntity * e) ;

	/** \brief return a copy of the behaviour
	 *
	 * @return a new Maxwell
	 */
	virtual Form * getCopy() const ;

	/** \brief return true*/
	virtual bool changed() const ;

	virtual void scale( double d )
	{
		param *= d ;
	}
	
	virtual Vector getForcesFromAppliedStress( Vector & data, Function & shape, const GaussPointArray & gp, std::valarray<Matrix> & Jinv, std::vector<Variable> & v) ;

	virtual Vector getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e, std::valarray<Matrix> & Jinv, std::vector<Variable> & v) ;
	
} ;

struct StandardLinearSolid : public LinearForm
{
	Matrix eta1 ;
	Matrix etaeq ;
	Matrix decay ;
	Matrix rig1 ;
	Matrix rig0 ;
	std::vector<Variable> v ;
	/** \brief Constructor
	*
	* @param rig Complete expression of the Cauchy-Green Strain Tensor
	* @param eta Complete expression of the viscosity Tensor
	*/
	StandardLinearSolid( const Matrix & rig0, const Matrix & rig1, const Matrix & eta ) ;

	virtual ~StandardLinearSolid() ;

	/** \brief Apply the law.
	 *
	 * The matrix is computed as: \f$ \nabla^T h_i K \nabla h_j + \dot{\nabla}^T h_i E \dot{\nabla} h_j\f$
	 * @param p_i first basis polynomial.
	 * @param p_j second basis polynomial.
	 * @param gp Gauss Points used for the quadrature
	 * @param Jinv Inverse Jacobian Matrices corresponding to the gauss points
	 * @param ret Matrix to store the result
	 * @param vm virtualMachine to use to compute the result
	 */
	virtual void apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const ;

	/** \brief is Element fractured
	 *
	 * @return false
	 */
	virtual bool fractured() const ;

	virtual ElementState * createElementState( IntegrableEntity * e) ;

	/** \brief return a copy of the behaviour
	 *
	 * @return a new Maxwell
	 */
	virtual Form * getCopy() const ;

	/** \brief return true*/
	virtual bool changed() const ;

	virtual void scale( double d )
	{
		param *= d ;
		etaeq *= d ;
		decay *= d ;
		rig0 *= d ;
		eta1 *= d ;
		rig1 *= d ;
	}
	
	virtual Vector getForcesFromAppliedStress( Vector & data, Function & shape, const GaussPointArray & gp, std::valarray<Matrix> & Jinv, std::vector<Variable> & v) ;

	virtual Vector getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e, std::valarray<Matrix> & Jinv, std::vector<Variable> & v) ;
	
} ;




	struct IterativeMaxwell : public LinearForm
	{
		double chartime ;
		std::vector<Variable> v ;
		
		std::vector<Vector> imposedStressAtGaussPoints ;
		
		double coeff_unext, coeff_uprev, coeff_aprev;
		
		IterativeMaxwell(const Matrix & rig, double eta) ;
		virtual ~IterativeMaxwell() ;

		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		virtual void step(double timestep, ElementState & currentState) ;
		virtual void updateElementState(double timestep, ElementState & currentState) const ;
		virtual void preProcess( double timeStep, ElementState & currentState ) ;
 		void preProcessAtGaussPoint(double timestep, ElementState & currentState, int g) ;	  
		
		virtual Form * getCopy() const ;
		virtual ElementState * createElementState( IntegrableEntity * e) ;

		virtual bool hasInducedForces() const { return true ; }

		virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
		virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
		virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
		
 		void setNumberOfGaussPoints(size_t n) ;
		virtual void getCoefficients(double timestep) ;
		void getInstantaneousCoefficients() ;
	} ;

} ;

#endif