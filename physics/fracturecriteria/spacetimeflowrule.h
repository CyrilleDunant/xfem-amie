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
#ifndef SPACE_TIME_FLOW_RULE_H__
#define SPACE_TIME_FLOW_RULE_H__

#include "fracturecriterion.h"
#include "maxstrain.h"
#include "../../mesher/delaunay_3d.h"
#include "../material_laws/material_laws.h"

namespace Amie {


/*PARSE SpaceTimeNonLocalDamageFlowRule FractureCriterion
	@string[file_name] // location of the file containing the linear interpolation of damage vs strain
*/
class SpaceTimeNonLocalDamageFlowRule : public MaximumStrain
{
public:
	LinearInterpolatedMaterialLaw * ruleData ;
	Function ruleFunction ;

	SpaceTimeNonLocalDamageFlowRule( const Function & f) ;
	SpaceTimeNonLocalDamageFlowRule( LinearInterpolatedMaterialLaw * data) ;
	SpaceTimeNonLocalDamageFlowRule( std::pair<Vector, Vector> data) ;
	SpaceTimeNonLocalDamageFlowRule( std::string file) ;

	virtual ~SpaceTimeNonLocalDamageFlowRule() { } ;

	virtual FractureCriterion * getCopy() const ;

	virtual double grade(ElementState &s)  ;
        virtual double gradeAtTime(ElementState &s, double t)  ;

} ;


}

#endif
