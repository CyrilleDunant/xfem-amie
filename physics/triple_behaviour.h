// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011

#ifndef __TRIPLE_BEHAVIOUR_H_
#define __TRIPLE_BEHAVIOUR_H_

#include "physics_base.h"

namespace Mu
{

/** \brief A Geometry determined triple behaviour.
 * Each Gauss point used is attributed a behaviour depending on whether it lies with respect with two concentric geometries determining three domains
 */
class TrimaterialInterface : public LinearForm
{
public:
	Geometry * inGeometry ;
	Geometry * outGeometry ;
	Form * inBehaviour ;
	Form * midBehaviour ;
	Form * outBehaviour ;
		
	Function xtransform ;
	Function ytransform ;
	
	/** \brief Constructor, set the Behaviour s and the delimiting Geometry s
	 * 
	 * @param in Geometry
	 * @param out Geometry
	 * @param inbehaviour 
	 * @param midbehaviour
	 * @param outbehaviour 
	 */
	TrimaterialInterface(Geometry * in,Geometry * out, Form * inbehaviour, Form * midbehaviour, Form * outbehaviour);
	
	virtual ~TrimaterialInterface();
	
	/** \brief Set the coordinate transformation functions to use to determin whether a point expressed in the local coordinates of an element lies in or out the Geometry
	 * 
	 * @param x x transformation
	 * @param y y transformation
	 */
	virtual void transform(const Function & x, const Function & y) ;
	
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
	
	virtual bool changed() const ;
	
	/** Check for fracture state
	 *
	 * @return true if the element is fractured
	 */
	virtual bool fractured() const ;
	
	/** get Copy of the behaviour
	 *
	 * @return pointer to the copy. Caller is responsible fior cleaning memory
	 */
	virtual Form * getCopy() const ;

	/** \brief Time-step the behaviours
	* This will step all the behaviours
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
	std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const; 

	/** \brief Acessor, return the highest fracture citerion*/
	virtual FractureCriterion * getFractureCriterion() const ;
	
	virtual void scale (double d) 
	{ 
		if(inBehaviour)
			inBehaviour->scale(d) ;
		if(midBehaviour)
			midBehaviour->scale(d) ;
		if(outBehaviour)
			outBehaviour->scale(d) ;
	}
	
} ;

} ;

#endif
