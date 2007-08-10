#ifndef __DUAL_BEHAVIOUR_H_
#define __DUAL_BEHAVIOUR_H_

#include "physics.h"

namespace Mu
{

class BimaterialInterface : public LinearForm
{
public:
	Matrix inTensor ;
	Geometry * inGeometry ;
	Form * inBehaviour ;
	Form * outBehaviour ;
		
	Function xtransform ;
	Function ytransform ;
	
	BimaterialInterface(Geometry * in, Form * inbehaviour, Form * outbehaviour);
	
	virtual ~BimaterialInterface();
	
	virtual void transform(const Function & x, const Function & y) ;
	
	virtual Matrix getTensor(const Point * p) const ;
	
	virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const ;
	
	virtual Matrix apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv) const ;
	
	virtual bool hasInducedForces() const ;
	
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
	
	virtual Vector getForces(const ElementState * s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const ;
	
} ;

} ;

#endif
