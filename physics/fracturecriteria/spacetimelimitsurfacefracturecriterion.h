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



typedef enum
{
    FRAME_CARTESIAN,
    FRAME_PRINCIPAL,
    FRAME_MODEL,
    FRAME_MATERIAL,
} ReferenceFrame ;


/*PARSE SpaceTimeLimitSurfaceFractureCriterion FractureCriterion
	@string[stress_measure] // function used to get a scalar measure of the stress
	@string[failure_surface] // function representing the stress-strain curve
	@string[requirements] // defines internal variables required by the failure surface
        @string<ReferenceFrame>[frame] FRAME_CARTESIAN // defines the frame in which the criterion is evaluated
        @string<bool>[positive] false // ignores negative components of the strains and stress tensors
*/
class SpaceTimeLimitSurfaceFractureCriterion : public FractureCriterion
{
protected:
	Function measure ;
	Function surface ;

        ReferenceFrame frame ;
        bool positive ;
	bool needStringVariable ;

	std::string surfaceYCoordinate ;
	std::string surfaceZCoordinate ;
	std::string surfaceTCoordinate ;
	std::string surfaceUCoordinate ;
	std::string surfaceVCoordinate ;
	std::string surfaceWCoordinate ;

	std::map<std::string, double> values ;



public:

        SpaceTimeLimitSurfaceFractureCriterion( Function m, Function s, ReferenceFrame frm = FRAME_PRINCIPAL, bool pos = false, std::string surface_y = std::string(), std::string surface_z = std::string(), std::string surface_t = std::string(), std::string surface_u = std::string(), std::string surface_v = std::string(), std::string surface_w = std::string()) ;
        SpaceTimeLimitSurfaceFractureCriterion( std::string m, std::string f, std::string requirements, ReferenceFrame frm = FRAME_PRINCIPAL, bool ps = false) ;

        virtual ~SpaceTimeLimitSurfaceFractureCriterion() {  }

	virtual FractureCriterion * getCopy() const ;

	virtual double grade(ElementState &s)  ;
        virtual double gradeAtTime(ElementState &s, double t)  ;

	virtual double getTensileLimit(const ElementState & s) const ;

} ;


}

#endif
