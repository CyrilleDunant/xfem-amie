//
// C++ Interface: radialstiffnessgradient
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __LIN_STIFF_GRADIENT
#define __LIN_STIFF_GRADIENT

#include "physics_base.h"

namespace Mu
{

/** \brief Linear elastic material with a stiffness gradient.*/
class LinearStiffnessGradient :public LinearForm
{
public:
	Matrix paramAlt ;
	Function s ;
	Point left ;
	Point right ;
	std::vector<Variable> v ;

	/** \brief Constructor, Build the two Cauchy-Green tensor for the interpolation
	 * 
	 * @param E_int First Young Modulus
	 * @param nu_int First Poisson Ratio
	 * @param E_ext Second Young's Modulus
	 * @param nu_ext Second Poisson Ratio
	 * @param l gradient starting point
	 * @param r gradient end point
	 */
	LinearStiffnessGradient(double E_int, double nu_int, double E_ext, double nu_ext, Point l, Point r) ;
	
	virtual ~LinearStiffnessGradient() ;
	
	/** \brief Set the coordinate transformation functions to use to determine the global location of the point at which the behaviour is evaluated
	 * 
	 * @param x x transformation
	 * @param y y transformation
	 */
	virtual void transform(const Function & x, const Function & y);
	
	/** \brief Return the stifness tensor at the point considered
		* The stiffness matrix is recomputed for each Gauss point using linear interpolation
		* @param p Point to check in local coordinates
	 */
	virtual Matrix getTensor(const Point & p) const ;
	
	/** \brief Apply the behaviour.
	* The stiffness matrix is recomputed for each Gauss point using linear interpolation 
	* @param p_i first shape function.
	* @param p_j second shape function.
	* @param gp Set of gauss points for numerical integration
	* @param Jinv inverse jacobian matrices at the gauss points
	* @param ret matrix in which to sore the results
	* @param vm pointer to the virtual machine dedicated for the computation
	*/
	virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
	
	/** \brief Check for fracture state
	 *
	 * @return true if the element is fractured
	 */
	virtual bool fractured() const ;
	
	/** \brief get Copy of the behaviour
	 *
	 * @return pointer to the copy. Caller is responsible fior cleaning memory
	 */
	virtual Form * getCopy() const ;
} ;

} ;

#endif
