//
// C++ Implementation: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "spacetimeflowrule.h"
#include "../damagemodels/damagemodel.h"
#include <fstream>

namespace Amie {

SpaceTimeNonLocalDamageFlowRule::SpaceTimeNonLocalDamageFlowRule( const Function & f, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain( 0., mirroring, delta_x, delta_y, delta_z), ruleData(nullptr), ruleFunction(f)
{

}

SpaceTimeNonLocalDamageFlowRule::SpaceTimeNonLocalDamageFlowRule( LinearInterpolatedMaterialLaw * data, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain( 0., mirroring, delta_x, delta_y, delta_z), ruleData(data), ruleFunction("0")
{

}

SpaceTimeNonLocalDamageFlowRule::SpaceTimeNonLocalDamageFlowRule( std::pair<Vector, Vector> data, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain( 0., mirroring, delta_x, delta_y, delta_z), ruleData(nullptr), ruleFunction("0")
{
   std::pair<std::string, std::string> tmp = std::make_pair( "strain", "damage") ;
   ruleData = new LinearInterpolatedMaterialLaw( tmp, data ) ;
}

SpaceTimeNonLocalDamageFlowRule::SpaceTimeNonLocalDamageFlowRule( std::string file, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain( 0., mirroring, delta_x, delta_y, delta_z), ruleData(nullptr), ruleFunction("0")
{
   std::pair<std::string, std::string> tmp = std::make_pair( "strain", "damage") ;
   ruleData = new LinearInterpolatedMaterialLaw( tmp, file ) ;
}

FractureCriterion * SpaceTimeNonLocalDamageFlowRule::getCopy() const 
{
    SpaceTimeNonLocalDamageFlowRule * ret = ( ruleData ? new SpaceTimeNonLocalDamageFlowRule( ruleData ) : new SpaceTimeNonLocalDamageFlowRule( ruleFunction ) ) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}

double SpaceTimeNonLocalDamageFlowRule::grade(ElementState &s)
{    
    double gradeBefore = gradeAtTime(s, -1) ;
    double gradeAfter = gradeAtTime(s, 1) ;
    scoreAtTimeStepEnd = gradeAfter ;

    if(gradeAfter < 0)
        return gradeAfter ;
    if(gradeBefore > 0)
    {
        return 1 ;
    }
    
    double upTime = 1 ;
    double downTime = -1 ;
    double testTime = 0 ;

    while(std::abs(upTime-downTime) > 1e-6)
    {
        double gradeTest = gradeAtTime(s, testTime) ;
        if(gradeTest < 0)
            downTime = testTime ;
        else if(gradeTest > 0)
            upTime = testTime ;
        else
            return testTime ;
        
        testTime = 0.5*(downTime+upTime) ;
    }
    return 1.-(testTime*.5+.5) ;
}


double SpaceTimeNonLocalDamageFlowRule::gradeAtTime(ElementState &s, double t)  
{
	Vector damage(1) ;
	Vector strain = getSmoothedField( PRINCIPAL_MECHANICAL_STRAIN_FIELD, s, t) ;
	s.getField( SCALAR_DAMAGE_FIELD, Point(0,0,0,t), damage, true ) ;

	if(s.getParent()->getBehaviour()->fractured() || damage[0] > 1.-POINT_TOLERANCE)
		return -1. ;

	double target = 0. ;
	if( ruleData )
		target = ruleData->get( strain.max() ) ;
	else
		target = VirtualMachine().eval( ruleFunction, strain.max() ) ;

	if(target < POINT_TOLERANCE)
		return -1. ;

	if( target > damage[0] )
		return std::min( 1., (target-damage[0])/target ) ;
	else
		return std::max(-1., -1.+(damage[0]-target)/damage[0]) ;
}


}

