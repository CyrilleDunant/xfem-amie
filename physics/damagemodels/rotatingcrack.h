//
// C++ Interface: isotropiclineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MU_ROTATING_CRACK_H
#define MU_ROTATING_CRACK_H

#include "damagemodel.h"

extern Mu::Matrix E;
namespace Mu {

/** \brief Rotating crack damage model. The stifness of an affected element is scaled by a factor between 1 and 0
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class RotatingCrack : public DamageModel
{
protected:
	double currentAngle ;
	double currentDamage ;
  
	bool inTension ;
	bool damaging ;
	
	double tdamage ;
	double cdamage ;
	
	double E ;
	double nu ;
	double factor ;
	
	std::vector< std::pair<double, double> > compressionAngles ;
	std::vector< std::pair<double, double> > tensionAngles ;
  std::vector<double> compressionweights ;
	std::vector<double> tensionweights ;
public:
	/** \brief Constructor. Set the number of degrees of freedom
	 * 
	 * @param numDof number of degrees of freedom
	 */
	RotatingCrack(double E, double nu) ;

	virtual ~RotatingCrack();
	virtual void scale(double s) { factor = s ;} ;

	/** \brief Increment the damage
	 * 
	 * @param s ElementState passed as a parameter
	 */
	virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;

	/** \brief compute the new stifness matrix after damage
	 * 
	 * \f$ K' = K(1-d) \f$
	 * @param m Matrix to modify
	 * @return the new Matrix
	 */
	virtual Matrix apply(const Matrix & m) const;
	
	virtual Matrix applyPrevious(const Matrix & m) const;

	/** \brief return true is the element concerned is fractured 
		*/
	virtual bool fractured() const  ;
	
	virtual void postProcess() ;
	
	virtual DamageModel * getCopy() const { return new RotatingCrack(E, nu) ;}
};

}

#endif
