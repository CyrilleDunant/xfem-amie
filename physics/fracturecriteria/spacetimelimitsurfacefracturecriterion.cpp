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
#include "spacetimelimitsurfacefracturecriterion.h"
#include "../damagemodels/damagemodel.h"
#include "../../utilities/parser.h"
#include <fstream>

namespace Amie {


SpaceTimeLimitSurfaceFractureCriterion::SpaceTimeLimitSurfaceFractureCriterion( Function m, Function s, StressMeasurementMethod mth, std::string y, std::string z, std::string t, std::string u, std::string v, std::string w, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z), measure(m), surface(s), method(mth), needStringVariable(false), surfaceYCoordinate(y), surfaceZCoordinate(z), surfaceTCoordinate(t), surfaceUCoordinate(u), surfaceVCoordinate(v), surfaceWCoordinate(w)
{
   if( surfaceYCoordinate.length() > 0 || surfaceZCoordinate.length() > 0 || surfaceTCoordinate.length() > 0 || surfaceUCoordinate.length() > 0 || surfaceVCoordinate.length() > 0 || surfaceWCoordinate.length() > 0 )
       needStringVariable = true ;
}

SpaceTimeLimitSurfaceFractureCriterion::SpaceTimeLimitSurfaceFractureCriterion( std::string m, std::string f, std::string reqs, StressMeasurementMethod mth, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion(mirroring, delta_x, delta_y, delta_z), method(mth), needStringVariable(false)
{
    std::map<std::string, char> coordMeasures ; 
    coordMeasures["stress_1"] = 'x' ;
    coordMeasures["stress_2"] = 'y' ;
    coordMeasures["stress_3"] = 'z' ;
    coordMeasures["stress_4"] = 't' ;
    coordMeasures["stress_5"] = 'u' ;
    coordMeasures["stress_6"] = 'v' ;
    measure = FunctionParser::getFunction( m, coordMeasures ) ;

    std::vector<std::string> req ;
    if(reqs.length() > 0)
    {
        size_t i = 0 ;
        size_t pos = reqs.find(',') ;
        while( pos != std::string::npos )
        {
            req.push_back( reqs.substr( i, pos-i ) ) ;
            i = ++pos ;
            pos = reqs.find(',', i) ;
        }
       req.push_back( reqs.substr(i, reqs.length() ) ) ;
    }
/*    for(size_t i = 0 ; i < req.size() ; i++)
	std::cout << "---" << req[i] << "---" << std::endl;*/

    std::map<std::string, char> coordLimit ;
    coordLimit["strain"] = 'x' ;
    if(req.size() > 0)
    {
        needStringVariable = true ;
        coordLimit[ req[0] ] = 'y' ;
        surfaceYCoordinate = req[0] ;
    }
    if(req.size() > 1)
    {
        needStringVariable = true ;
        coordLimit[ req[1] ] = 'z' ;
        surfaceZCoordinate = req[1] ;
    }
    if(req.size() > 2)
    {
        needStringVariable = true ;
        coordLimit[ req[2] ] = 't' ;
        surfaceTCoordinate = req[2] ;
    }
    if(req.size() > 3)
    {
        needStringVariable = true ;
        coordLimit[ req[3] ] = 'u' ;
        surfaceUCoordinate = req[3] ;
    }
    if(req.size() > 4)
    {
        needStringVariable = true ;
        coordLimit[ req[4] ] = 'v' ;
        surfaceVCoordinate = req[4] ;
    }
    if(req.size() > 5)
    {
        needStringVariable = true ;
        coordLimit[ req[5] ] = 'w' ;
        surfaceWCoordinate = req[5] ;
    }
    surface = FunctionParser::getFunction( f, coordLimit ) ;

}

FractureCriterion * SpaceTimeLimitSurfaceFractureCriterion::getCopy() const 
{
    SpaceTimeLimitSurfaceFractureCriterion * ret = new SpaceTimeLimitSurfaceFractureCriterion( measure, surface, method, surfaceYCoordinate, surfaceZCoordinate, surfaceTCoordinate, surfaceUCoordinate, surfaceVCoordinate, surfaceWCoordinate ) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}

double SpaceTimeLimitSurfaceFractureCriterion::grade(ElementState &s)
{    
    double gradeBefore = gradeAtTime(s, -1) ;
    double gradeAfter = gradeAtTime(s, 1) ;
    scoreAtTimeStepEnd = gradeAfter ;

    if(gradeAfter < 0)
        return gradeAfter ;
    if(gradeBefore > -POINT_TOLERANCE)
        return 1 ;
    
    double upTime = 1 ;
    double downTime = -1 ;
    double testTime = 2.*(-gradeBefore)/(gradeAfter-gradeBefore)-1. ;
    double gradeDown = gradeBefore ;
    double gradeUp = gradeAfter ;
    double gradeTest = 0 ;

    while(std::abs(upTime-downTime) > 1e-4 )
    {
        gradeTest = gradeAtTime(s, testTime) ;
        if(gradeTest < 0)
        {
            downTime = testTime ;
            gradeDown = gradeTest ;
        }
        else if(gradeTest > 0)
        {
            upTime = testTime ;
            gradeUp = gradeTest ;
        }
        else
            return 1.-(testTime*.5+.5) ;
        
	if(gradeDown > POINT_TOLERANCE)
	        testTime = downTime + (upTime-downTime)*(-gradeDown)/(gradeUp-gradeDown) ;
	else
		testTime = (downTime + upTime)*0.5 ;
	if(std::abs(testTime - upTime) < POINT_TOLERANCE && std::abs(gradeTest) < POINT_TOLERANCE) 
		return 1.-(testTime*.5+.5) ;
    }

    return 1.-(testTime*.5+.5) ;
}

double getVal( Vector & v, size_t iter)
{
    if(v.size() > iter) { return v[iter] ; }
    return 0 ;
}

double SpaceTimeLimitSurfaceFractureCriterion::gradeAtTime(ElementState &s, double t)  
{
    if( s.getParent()->getBehaviour()->fractured() )
        return -1 ;

    FieldType sigma = REAL_STRESS_FIELD ;
    FieldType epsilon = MECHANICAL_STRAIN_FIELD ;
    if( method == PRINCIPAL || method == PRINCIPAL_POSITIVE || method == PRINCIPAL_NEGATIVE )
    {
        sigma = PRINCIPAL_REAL_STRESS_FIELD ;
        epsilon = PRINCIPAL_MECHANICAL_STRAIN_FIELD ;
    }
    std::pair<Vector, Vector> currentState = getSmoothedFields( sigma, epsilon, s, t ) ;
    switch( method )
    {
        case ALL:
        case PRINCIPAL:
            break;
        case ALL_POSITIVE:
        case PRINCIPAL_POSITIVE:
            for(size_t i = 0 ; i < currentState.first.size() ; i++)
                currentState.first[i] = std::max(0., currentState.first[i] ) ;
            for(size_t i = 0 ; i < currentState.second.size() ; i++)
                currentState.second[i] = std::max(0., currentState.second[i] ) ;
            break;
        case ALL_NEGATIVE:
        case PRINCIPAL_NEGATIVE:
            for(size_t i = 0 ; i < currentState.first.size() ; i++)
                currentState.first[i] = std::max(0., -currentState.first[i] ) ;
            for(size_t i = 0 ; i < currentState.second.size() ; i++)
                currentState.second[i] = std::max(0., -currentState.second[i] ) ;
            break;
    }

    VirtualMachine vm ;
    double stress = vm.eval( measure, getVal( currentState.first, 0), getVal( currentState.first, 1), getVal( currentState.first, 2), getVal( currentState.first, 3), getVal( currentState.first, 4), getVal( currentState.first, 5)) ;
    double strain = vm.eval( measure, getVal( currentState.second, 0), getVal( currentState.second, 1), getVal( currentState.second, 2), getVal( currentState.second, 3), getVal( currentState.second, 4), getVal( currentState.second, 5)) ;

    double ycoor = 0 ;
    double zcoor = 0 ;
    double tcoor = 0 ;
    double ucoor = 0 ;
    double vcoor = 0 ;
    double wcoor = 0 ;

    if( needStringVariable )
    {
        if( surfaceYCoordinate.length() > 0 )
            ycoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( surfaceYCoordinate, values ) ;
        if( surfaceZCoordinate.length() > 0 )
            zcoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( surfaceZCoordinate, values ) ;
        if( surfaceTCoordinate.length() > 0 )
            tcoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( surfaceTCoordinate, values ) ;
        if( surfaceUCoordinate.length() > 0 )
            ucoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( surfaceUCoordinate, values ) ;
        if( surfaceVCoordinate.length() > 0 )
            vcoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( surfaceVCoordinate, values ) ;
        if( surfaceWCoordinate.length() > 0 )
            wcoor = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get( surfaceWCoordinate, values ) ;
    }


    double limit = std::max(0., vm.eval( surface, strain, ycoor, zcoor, tcoor, ucoor, vcoor, wcoor )) ;

/*    std::cout << strain << "\t" << stress << "\t" << limit << std::endl ;*/

    if( stress > limit )
        return std::min( 1., 1.-limit/stress ) ;

    return std::max( -1., -1.+ stress/limit ) ;
}

double SpaceTimeLimitSurfaceFractureCriterion::getTensileLimit(const ElementState & s) const 
{
    VirtualMachine vm ;

    double ycoor = 0 ;
    double zcoor = 0 ;
    double tcoor = 0 ;
    double ucoor = 0 ;
    double vcoor = 0 ;
    double wcoor = 0 ;

    if( needStringVariable )
    {
        std::map<std::string, double> copied = dynamic_cast<const GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).getVariables() ;
        if( surfaceYCoordinate.length() > 0 )
            ycoor = copied[surfaceYCoordinate] ;
        if( surfaceZCoordinate.length() > 0 )
            zcoor = copied[surfaceZCoordinate] ;
        if( surfaceTCoordinate.length() > 0 )
            tcoor = copied[surfaceTCoordinate] ;
        if( surfaceUCoordinate.length() > 0 )
            ucoor = copied[surfaceUCoordinate] ;
        if( surfaceVCoordinate.length() > 0 )
            vcoor = copied[surfaceVCoordinate] ;
        if( surfaceWCoordinate.length() > 0 )
            wcoor = copied[surfaceWCoordinate] ;
    }

    return std::max(0., vm.eval( 0., ycoor, zcoor, tcoor, ucoor, vcoor, wcoor )) ;

}


/*std::string SpaceTimeLimitSurfaceFractureCriterion::getStressMeasurementFunction( StressMeasurementFunction f, SpaceDimensionality dim) 
{
    switch(f)
    {
        case MAXIMUM_PRINCIPAL_STRESS:
        {
            if(dim == SPACE_THREE_DIMENSIONAL)
                return "max stress_1 ( max stress_2 stress_3 )" ;
            return "max stress_1 stress_2" ;
        }
        case FIRST_STRESS_INVARIANT:
        {
            if(dim == SPACE_THREE_DIMENSIONAL )
                return "stress_1 + stress_2 + stress_3" ;
            return "stress_1 + stress_2" ;
        }
        case SECOND_STRESS_INVARIANT:
        {
            if(dim == SPACE_THREE_DIMENSIONAL )
                return "sqrt ( stress_1 * stress_2 + stress_2 * stress_3 + stress_3 * stress_1 - stress_4 * stress_4 - stress_5 * stress_5 - stress_6 * stress_6 )" ;
            return "sqrt ( stress_1 * stress_2 - stress_3 * stress_3 )" ;
        }
        case THIRD_STRESS_INVARIANT:
        {
            if(dim == SPACE_THREE_DIMENSIONAL )
                return "( stress_1 * stress_2 * stress_3 + 2 * stress_4 * stress_5 * stress_6 - stress_6 * stress_6 * stress_3 - stress_4 * stress_4 * stress_1 - stress_5 * stress_5 * stress_2 ) ^ 0.3333333333333333" ;
            return "0" ;
        }
        case SECOND_DEVIATORIC_STRESS_INVARIANT:
        {
            if(dim == SPACE_THREE_DIMENSIONAL )
                return "sqrt ( ( ( stress_1 - stress_2 ) * ( stress_1 - stress_2 ) + ( stress_2 - stress_3 ) * ( stress_2 - stress_3 ) + ( stress_3 - stress_1 ) * ( stress_3 - stress_1 ) ) / 6 )" ;
            return "sqrt ( ( ( stress_1 - stress_2 ) * ( stress_1 - stress_2 ) ) / 6 )" ;
        }
        case THIRD_DEVIATORIC_STRESS_INVARIANT:
        {
            if(dim == SPACE_THREE_DIMENSIONAL )
                std::cerr << "warning, space measurement function not implemented!" << std::endl ;
                return "0" ;
            return "0" ;
        }
        case MAXIMUM_SHEAR_STRESS:
        {
            if(dim == SPACE_THREE_DIMENSIONAL )
                return "0.5 * max ( abs ( stress_1 - stress_ 2) ) ( max ( abs ( stress_2 - stress_3 ) ) ( abs ( stress_3 - stress_1 ) ) )" ;
            return "0.5 * max ( abs ( stress_1 - stress_ 2) ) ( max ( abs ( stress_2 ) ) ( abs ( stress_1 ) ) )" ;
        }
        case VON_MISES_STRESS:
        {
            if(dim == SPACE_THREE_DIMENSIONAL )
                return "sqrt ( ( stress_1 - stress_2 ) * ( stress_1 - stress_2 ) + ( stress_2 - stress_3 ) * ( stress_2 - stress_3 ) + ( stress_3 - stress_1 ) * ( stress_3 - stress_1 ) )" ;
            return "sqrt ( ( stress_1 - stress_2 ) * ( stress_1 - stress_2 ) + stress_1 * stress_1 + stress_2 * stress_2 ) " ;
        }
    }
    if(dim == SPACE_THREE_DIMENSIONAL)
        return "max stress_1 ( max stress_2 stress_3 )" ;
    return "max stress_1 stress_2" ;
}*/



}

