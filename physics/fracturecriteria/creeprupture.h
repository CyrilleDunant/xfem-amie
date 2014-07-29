//
// C++ Interface: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef CREEPRUPTURE_H
#define CREEPRUPTURE_H

#include "fracturecriterion.h"

namespace Amie {

/** \brief Mohr - Coulomb fracture criterion
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The Mohr-Coulomb is met when the maximum principal stresses is below or above the prescribed limits 
	
*/

class CreepRupture : public FractureCriterion
{
public:
	double maxStress ;
	double limStress ;
	double limStrain ;
	bool metInCompression  ;
	bool metInTension  ;
	
	virtual bool directionInTension(size_t direction) {return metInCompression ;}
	virtual bool directionInCompression(size_t direction) {return metInTension ;}

	CreepRupture(double maxStress, double limStress, double limStrain, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);

	virtual ~CreepRupture();

	virtual FractureCriterion * getCopy() const;

	virtual double grade(ElementState &s)  ;

	virtual Material toMaterial() ;
		
	virtual double getTensileLimit(const ElementState & s) const {return maxStress ;};

};

}

#endif
