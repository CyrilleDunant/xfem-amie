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

#ifndef __STIFFNESS_H_
#define __STIFFNESS_H_

#include "physics_base.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/damagemodel.h"
#include "homogenization/homogenization_base.h"
#include "fracturecriteria/nonlocalvonmises.h"

namespace Mu
{

	/** \brief A linear Elastic Law
	* The field param is the Cauchy-Green Strain Tensor
	*/
	struct Stiffness : public LinearForm
	{
		std::vector<Variable> v ;
		/** \brief Constructor
		* 
		* @param rig Complete expression of the Cauchy-Green Strain Tensor
		*/
		Stiffness(const Matrix & rig) ;
		
		virtual ~Stiffness() ;

		virtual XMLTree * toXML() { return new XMLTree("stiffness",param) ; } ;
		
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
		
		virtual void stepBack() { };
		
		/** \brief Return a copy of the behaviour*/
		virtual Form * getCopy() const ;

	} ;

	
	/** \brief A linear pseudo-plastic law.
	* The behaviour is updated when step is called and the secant stifness is computed as a function  of the strain
	*/
	class PseudoPlastic : public LinearForm
	{
	protected:
		void fixLastDamage() ;
		NonLocalVonMises  * vm ;
		bool initialised ;
		
	public:
		std::vector<DelaunayTriangle *> cache ;
		std::vector<Variable> v ;
		double alpha ;
		double lastDamage ;
		bool change ;
		bool frac ;
		bool fixedfrac ;
		double radius ;
		double limitStrain ;
		
		/** \brief Constructor
		* 
		* @param rig Complete expression of the Cauchy-Green Strain Tensor
		*/
		PseudoPlastic(const Matrix & rig, double limitStrain, double radius) ;
		
		
		
		virtual ~PseudoPlastic() ;

		virtual XMLTree * toXML() 
		{ 
			return new XMLTree("pseudoplastic",param) ; 
			
		} ;
		
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
		
		/** \brief update behaviour
		*
		* @param timestep elapsed time
		* @param currentState state of the element
		* 
		*/
		virtual void step(double timestep, ElementState & currentState) ;
		
		virtual bool changed() const ;
		virtual Matrix getTensor(const Point & p) const ;
		virtual Matrix getPreviousTensor(const Point & p) const ;
		virtual FractureCriterion * getFractureCriterion() const ;
		
		/** \brief Return false.*/
		virtual bool fractured() const ;
		
		/** \brief Return a copy of the behaviour*/
		virtual Form * getCopy() const ;
		
		
	} ;
	
} ;

#endif
