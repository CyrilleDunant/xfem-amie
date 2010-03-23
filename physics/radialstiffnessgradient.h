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

#ifndef __RADIAL_STIFF_GRADIENT
#define __RADIAL_STIFF_GRADIENT

#include "physics_base.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/lineardamage.h"

namespace Mu
{

/** \brief  Radial stiffness gradient: the mechanical property is scaled according to the position between two radii. */
class RadialStiffnessGradient :public LinearForm
{
public:
	Matrix paramAlt ;
	Function r ;
	double r_ext ;
	double r_int ;
	double dr ;
	
	double previousDamage ;
	FractureCriterion * criterion ;
	bool frac ;
	bool change ; 
	double sigmaRupt ;
	double init ;
	double damage ;
	Point centre ;
	LinearDamage dfunc ;

	/** \brief  Radial stiffness gradient: the mechanical property is scaled according to the position between two radii.
	*
	* No FractureCriterion is set at initialisation: the material is "unbreakable"
	* @param E_int Interior Young's Modulus
	* @param nu_int Interior Poisson ratio
	* @param rint Interior radius
	* @param E_ext Exterior Young's Modulus
	* @param nu_ext Exterior Poisson ratio
	* @param rext exterior radius
	* @param c Center
*/ 
	RadialStiffnessGradient(double E_int, double nu_int, double rint, double E_ext, double nu_ext, double rext, Point c) ;
	
	virtual ~RadialStiffnessGradient() ;

		/** \brief Set the coordinate transformation functions to use to determine the global location of the point at which the behaviour is evaluated
	 * 
	 * @param x x transformation
	 * @param y y transformation
	 */
	virtual void transform(const Function & x, const Function & y);
	
	/** \brief Apply the behaviour.
	*
		* The stiffness matrix is recomputed for each Gauss point using linear interpolation
		* @param p_i first basis polynomial.
		* @param p_j second basis polynomial.
		* @param e IntegrableEntity on which to perform the integration.
	 */
	virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const ;
	
	/** \brief Return the stifness tensor at the point considered
	*
		* The stiffness matrix is recomputed for each Gauss point using linear interpolation
		* @param p Point to check in local coordinates
	 */
	virtual Matrix getTensor(const Point & p) const ;
	
	/** \brief Apply the behaviour.
	*
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
	
	virtual void getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector &v) const ;
	
	/** \brief Check for fracture
	* @param timestep elapsed time
	* @param currentState state of the element
	* 
	* if the Criterion is true, se fractured state to true
	*/
	virtual void step(double timestep, ElementState & currentState) ;

	virtual void artificialDamageStep(double d) ;

	/** \brief return true if the damage state changed during the last step */
	virtual bool changed() const ;

/** \brief Unwind a step in the behaviour history
 *
 * The damage state will return to its previous state
*/
	virtual void stepBack() ;

/** \brief Set a FractureCriterion 
 *
 * The behaviour will delete the Criterion when it is itself destroyed
 * @param crit New Criterion to set
*/
	void setFractureCriterion(FractureCriterion * crit) ;
} ;

} ;

#endif
