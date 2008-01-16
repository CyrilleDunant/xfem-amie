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

namespace Mu
{

class RadialStiffnessGradient :public LinearForm
{
public:
	Matrix paramAlt ;
	Function r ;
	double r_ext ;
	double r_int ;
	double dr ;
	
	Point centre ;

	RadialStiffnessGradient(double E_int, double nu_int, double rint, double E_ext, double nu_ext, double rext, Point c) ;
	
	virtual ~RadialStiffnessGradient() ;
	
	virtual void transform(const Function & x, const Function & y);
	
	virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const ;
	
	virtual Matrix getTensor(const Point & p) const ;
	
	virtual Matrix apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
	
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
	
	virtual Vector getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
	
} ;

} ;

#endif
