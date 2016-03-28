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
#ifndef SPACE_TIME_LIMITSURFACE_FRACTURE_CRITERION_H__
#define SPACE_TIME_LIMITSURFACE_FRACTURE_CRITERION_H__

#include "fracturecriterion.h"

namespace Amie {

typedef enum {
    ALL,
    ALL_POSITIVE,
    ALL_NEGATIVE,
    PRINCIPAL,
    PRINCIPAL_POSITIVE,
    PRINCIPAL_NEGATIVE,
} StressMeasurementMethod ;




class SpaceTimeLimitSurfaceFractureCriterion : public FractureCriterion
{
protected:
	Function measure ;
	Function surface ;

	StressMeasurementMethod method ;
	bool needStringVariable ;

	std::string surfaceYCoordinate ;
	std::string surfaceZCoordinate ;
	std::string surfaceTCoordinate ;
	std::string surfaceUCoordinate ;
	std::string surfaceVCoordinate ;
	std::string surfaceWCoordinate ;

	std::map<std::string, double> values ;


public:

	SpaceTimeLimitSurfaceFractureCriterion( Function m, Function s, StressMeasurementMethod mth = ALL, std::string surface_y = std::string(), std::string surface_z = std::string(), std::string surface_t = std::string(), std::string surface_u = std::string(), std::string surface_v = std::string(), std::string surface_w = std::string(), MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0 ) ;
	SpaceTimeLimitSurfaceFractureCriterion( std::string m, std::string f, std::string requirements, StressMeasurementMethod mth = ALL, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0 ) ;

	virtual ~SpaceTimeLimitSurfaceFractureCriterion() {  } ;

	virtual FractureCriterion * getCopy() const ;

	virtual double grade(ElementState &s)  ;
        virtual double gradeAtTime(ElementState &s, double t)  ;

	virtual double getTensileLimit(const ElementState & s) const ;

} ;


}

#endif
