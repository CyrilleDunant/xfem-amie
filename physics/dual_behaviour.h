#ifndef __DUAL_BEHAVIOUR_H_
#define __DUAL_BEHAVIOUR_H_

#include "physics_base.h"

namespace Mu
{

/** \brief A Geometry determined dual behaviour.
 * Each Gauss point used is attributed a behaviour depending on whether it lies in or out a given Geometry
 */
class BimaterialInterface : public LinearForm
{
public:
	Geometry * inGeometry ;
	Form * inBehaviour ;
	Form * outBehaviour ;
		
	Function xtransform ;
	Function ytransform ;
	Function ztransform ;
	
	/** \brief Constructor, set the Behaviour s and the delimiting Geometry
	 * 
	 * @param in Geometry
	 * @param inbehaviour 
	 * @param outbehaviour 
	 */
	BimaterialInterface(Geometry * in, Form * inbehaviour, Form * outbehaviour);
	
	virtual ~BimaterialInterface();
	
	/** \brief Set the coordinate transformation functions to use to determin whether a point expressed in the local coordinates of an element lies in or out the Geometry
	 * 
	 * @param x x transformation
	 * @param y y transformation
	 */
	virtual void transform(const Function & x, const Function & y) ;

	/** \brief Set the coordinate transformation functions to use to determin whether a point expressed in the local coordinates of an element lies in or out the Geometry
	 * 
	 * @param x x transformation
	 * @param y y transformation
	 * @param z z transformation
	 */
	virtual void transform(const Function & x, const Function & y, const Function & z) ;
	
	/** \brief return the linear factor of the behaviour corresponding to the position given
	 * 
	 * @param p Point at which to compute the linear parameter
	 * @return 
	 */
	virtual Matrix getTensor(const Point & p) const ;

	/** \brief Return the imposed stress at the point considered
	 * 
	 * @param p Point
	 * @return stress Vector 
	 */
	virtual Vector getImposedStress(const Point & p) const ;
	
/** \brief Apply the behaviour
	* This overloaded apply() is more efficient and is designed to minimise allocating and dealocating memory.
	* The result of the computation depends on the location of the Gauss points
	* 
	* @param p_i first shape function.
	* @param p_j second shape function.
	* @param gp Set of gauss points for numerical integration
	* @param Jinv inverse jacobian matrices at the gauss points
	* @param ret matrix in which to sore the results
	* @param vm pointer to the virtual machine dedicated for the computation
	*/
	virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
	
	/** \brief Check if this behaviour induces internal forces
	 * 
	 * @return true if any of the two behaviour has induces forces
	 */
	virtual bool hasInducedForces() const ;
	
	/** \brief Check for fracture state
	 *
	 * @return true if the element is fractured
	 */
	virtual bool fractured() const ;
	
	/** \brief get a copy of the behaviour
	 * This will create a new Bimateral behaviour, with a copy of both the members.
	 * @return pointer to the copy. Caller is responsible fior cleaning memory
	 */
	virtual Form * getCopy() const ;

	/** \brief Time-step the behaviour
	* This will step both be behaviours
	* @param timestep delta-time of the step.
	* @param currentState State of the element in which the behaviour is time-stepped
	*/
	virtual void step(double timestep, ElementState & currentState) ;

	virtual void artificialDamageStep(double d) ;
	
	/** \brief Return the vector of induced forces if any of the behaviours induces internal forces. Return an empty vecor otherwise
	 * 
	 * @param s ElementState to consider
	 * @param p_i shape function to consider
	 * @param gp Gauss points to use for the computation of the quadrature
	 * @param Jinv Inverse Jacobian matrix to use at the Gauss Points
	 * @param v Vector in which to store the result
	 */
	virtual void getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector &v) const ;
	
} ;

} ;

#endif
