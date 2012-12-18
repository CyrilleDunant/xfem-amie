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



	/** \brief generic class for iterative Maxwell behaviour. This class can be inherited to provide alternative iteration schemes */
	struct IterativeMaxwell : public LinearForm
	{
		/** \brief characteristic time of the branch */
		double chartime ;
		std::vector<Variable> v ;
		/** \brief list of viscoelastic imposed stresses at all Gauss points */
		std::vector<Vector> imposedStressAtGaussPoints ;
		/** \brief coefficients relating a_{n+1} to u_{n+1} u_n and a_n. These coefficients are updated at each time step to account for the state of the system */
		double coeff_unext, coeff_uprev, coeff_aprev;
		
		/** \brief Constructor 
		 * @param rig stiffness matrix of the spring of the branch
		 * @param eta characteristic time of the dashpot of the branch
		 */
		IterativeMaxwell(const Matrix & rig, double eta) ;
		virtual ~IterativeMaxwell() ;

		/** \brief construct the elementary matrix associated to the nodes i and j in an element
		 * @param p_i shape function
		 * @param p_j shape function
		 * @param gp Gauss points of the element
		 * @param Jinv list of inverse jacobian matrices of the element at each Gauss point
		 * @param ret storage matrix
		 * @param vm virtual machine
		 */
		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		/** \brief updates the behaviour (does nothing) */
		virtual void step(double timestep, ElementState & currentState) ;
		/** \brief updates the element state according to the current state
		 * @param timestep duration of the last time step
		 * @param currentState state of the current element after the last time step
		 */
		virtual void updateElementState(double timestep, ElementState & currentState) const ;
		/** \brief modifies the rigidity tensor and the viscoelastic stresses to account for the state of the current element 
		 * @param timestep the duration of the upcoming time step
		 * @param currentState the element state of the current element after the last time step
		 */
		virtual void preProcess( double timeStep, ElementState & currentState ) ;
		/** \brief modifies the imposed stress at the current Gausse point
		 * @param timestep the duration of the upcoming time step
		 * @param currentState the element state of the current element after the last time step
		 * @param g index of the Gauss point
		 */
 		void preProcessAtGaussPoint(double timestep, ElementState & currentState, int g) ;	  
		
		/** \brief instantiate a hard copy of the current behaviour */
		virtual Form * getCopy() const ;

		/** \brief instantiate an element state containing internal variables corresponding to the strain and the anelastic strain at the previous time step 
		 * @param e element		 
		 */
		virtual ElementState * createElementState( IntegrableEntity * e) ;

		/** \brief flag to force application of boundary conditions */
		virtual bool hasInducedForces() const { return true ; }

		/** \brief returns the imposed stress at a given point
		 * @param p point
		 * @param e element for the discretization
		 * @param g index of Gauss point (if p is a Gauss point)
		 */
		virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
		
		/** \brief returns the imposed strain at a given point
		 * @param p point
		 * @param e element for the discretization
		 * @param g index of Gauss point (if p is a Gauss point)
		 */
		virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
		
		/** \brief transforms the imposed stress in an element into boundary conditions that are applied to the right-hand side of the equation
		 * @param s state of the element before the time step
		 * @param id node at which the boundary conditions are calculated
		 * @param p_i shape function associated to the node
		 * @param gp Gauss points of the element
		 * @param Jinv inverse jacoblian matrices calculated at each Gauss points
		 */
		virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
		
		/** \brief adapts the number of Gauss points
		 * @param n the new number of Gauss points
		 */
		void setNumberOfGaussPoints(size_t n) ;
		
		/** \brief creates the coefficients for the iteration. The coefficients are given in this format:
		 * a_{n+1} = coeff_unext u_{n+1} + coeff_uprev u_n + coeff_aprev a_n
		 * they are used to modify the stiffness tensor of the system and the imposed stress at each Gauss points
		 * You can create various iterative procedure by inheriting this structure and overload this method
		 * @param timestep duration of the upcoming time step
		 */
		virtual void getCoefficients(double timestep) ;

		/** \brief creates the coefficients for the iteration for instantaneous time step : a_{n+1} = a_n */
		void getInstantaneousCoefficients() ;
	} ;

	
	
	struct GeneralizedIterativeMaxwell : public LinearForm
	{
		std::vector<Variable> v ;
		std::vector<Vector> imposedStressAtGaussPoints ;
		std::vector<IterativeMaxwell *> branches ;
		Matrix r0 ;
		
		GeneralizedIterativeMaxwell(const Matrix & rig, const std::vector<Matrix> & param, const std::vector<double> & chartime) ;
		virtual ~GeneralizedIterativeMaxwell() ;

		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
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
		void syncNumberOfGaussPoints(ElementState & state) ;
		virtual void getCoefficients(double timestep) ;
		void getInstantaneousCoefficients() ;
	} ;
	

} ;

#endif