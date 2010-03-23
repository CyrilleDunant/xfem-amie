//
// C++ Interface: RadialDistributedStiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __RADIAL_DISTRIBUTED_RadialDistributedStiffness_H_
#define __RADIAL_DISTRIBUTED_RadialDistributedStiffness_H_

#include "physics_base.h"
#include "../features/inclusion.h"
#include "../utilities/xml.h"

namespace Mu
{

	/** \brief A linear Elastic Law
	* The field param is the Cauchy-Green Strain Tensor
	*/
	struct RadialDistributedStiffness : public LinearForm
	{
		std::vector<std::pair<double, Matrix> > stiff ;
		double angle ;

		std::vector<Variable> v ;
		/** \brief Constructor
		* 
		* @param rig Complete expression of the Cauchy-Green Strain Tensor
		*/
		RadialDistributedStiffness(std::vector<std::pair<double, Matrix> > rig) ;

		RadialDistributedStiffness(XMLTree * xml) ;
		
		virtual ~RadialDistributedStiffness() ;

		virtual XMLTree * toXML() ;

		void setAngle(double a) ;
		
		/** \brief Apply the law.
		* 
		* @param p_i first basis polynomial.
		* @param p_j second basis polynomial.
		* @return symbolic matrix resulting of \f$ \nabla H^T K \nabla H \f$.
		*/
		virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const; 
		
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
		
		/** \brief Return a 0-length Vector*/
		virtual void getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector &v) const ;
		
	} ;

	class RadialInclusion : public Inclusion
	{
		public:
			RadialInclusion(Feature * f, double r, Point c) ;
			RadialInclusion(Feature * f, double r, double x, double y) ;
			RadialInclusion(double r, Point c) ;
			RadialInclusion(double r, double x, double y) ;
			RadialInclusion(XMLTree * xml) ;

			virtual Form * getBehaviour( const Point & p ) ;		
	} ;	

} ;

#endif
