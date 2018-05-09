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
#ifndef SPACE_TIME_MULTISURFACE_FRACTURE_CRITERION_H__
#define SPACE_TIME_MULTISURFACE_FRACTURE_CRITERION_H__

#include "fracturecriterion.h"

namespace Amie {


/*PARSE SpaceTimeMultiSurfaceFractureCriterion FractureCriterion */
class SpaceTimeMultiSurfaceFractureCriterion : public FractureCriterion
{
protected:
	std::vector<FractureCriterion *> surfaces ;
	std::vector<double> instants ;

public:

	SpaceTimeMultiSurfaceFractureCriterion( ) : FractureCriterion() { }
	void add( FractureCriterion * crit ) { surfaces.push_back(crit) ; instants.push_back(-1) ; }

	virtual ~SpaceTimeMultiSurfaceFractureCriterion() {  } ;

	virtual FractureCriterion * getCopy() const ;

	virtual double grade(ElementState &s)  ;
        virtual double gradeAtTime(ElementState &s, double t)  ;

	virtual void initialiseCache( ElementState& s ) ;
	virtual void updateCache( ElementState & s);

	virtual bool directionMet(size_t direction, double t = 0) ;


} ;


}

#endif
