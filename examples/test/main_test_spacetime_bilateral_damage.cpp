// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/viscoelasticity.h"
#include "../../physics/viscoelasticity_and_fracture.h"
#include "../../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
#include "../../physics/fracturecriteria/spacetimemultisurfacefracturecriterion.h"
#include "../../physics/fracturecriteria/spacetimelimitsurfacefracturecriterion.h"
#include "../../physics/damagemodels/spacetimefiberbasedbilineardamage.h"
#include "../../physics/materials/paste_behaviour.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../features/sample.h"
#include "../../utilities/parser/command_line_parser.h"

#include <fstream>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

int main(int argc, char *argv[])
{
	CommandLineParser parser("Test a space-time behaviour in bulk and in shear") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::string outdir = parser.getString("--output-directory") ;


        Sample rect(nullptr, 0.01,0.01,0,0) ;
	SpaceTimeMultiSurfaceFractureCriterion * crit = new SpaceTimeMultiSurfaceFractureCriterion() ;
	SpaceTimeLimitSurfaceFractureCriterion * bulk = new SpaceTimeLimitSurfaceFractureCriterion( "( stress_1 + stress_2 ) / 2", "max ( bulk_strength ) ( bulk_strength * ( strain * bulk_modulus / bulk_strength ) )", "bulk_strength,bulk_modulus", PRINCIPAL ) ;
	SpaceTimeLimitSurfaceFractureCriterion * shear = new SpaceTimeLimitSurfaceFractureCriterion( "( sqrt ( ( stress_1 - stress_2 ) ^ 2 + stress_1 ^ 2 + stress_2 ^ 2 ) ) / 2", "shear_strength", "shear_strength", PRINCIPAL ) ;
	crit->add( bulk ) ;
	crit->add( shear ) ;
	SpaceTimeFiberBasedBilateralLinearDamage * dam = new SpaceTimeFiberBasedBilateralLinearDamage( 0.01, 1e-6, 0.99999 ) ;

	rect.setBehaviour(new LogarithmicCreepWithExternalParameters( "bulk_modulus = 10e9, shear_modulus = 5e9, bulk_strength = 1e6, shear_strength = 2e6", crit, dam ) ) ;

	FeatureTree f(&rect) ;
	f.setSamplingNumber(0) ;
        f.setDeltaTime(0.001) ;
        f.setMinDeltaTime(1e-9) ;

	f.step() ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
	BoundingBoxDefinedBoundaryCondition * top = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, 0) ;
	BoundingBoxDefinedBoundaryCondition * right = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_XI, TOP_RIGHT_AFTER, 0) ;
	f.addBoundaryCondition( top ) ;
//	f.addBoundaryCondition( right ) ;

	std::ofstream out ;
	if(renew)
		out.open(outdir+"/test_spacetime_bilateral_damage_base", std::ios::out) ;
	else
		out.open(outdir+"/test_spacetime_bilateral_damage_current", std::ios::out) ;

	for(size_t i = 0 ; i < 40 ; i++)
	{
		top->setData( 0.0000009 + 0.0000001*i ) ;
		right->setData( -0.5*(0.0000009 + 0.0000001*i) ) ;
		f.step() ;
		Vector stress = f.getAverageField( PRINCIPAL_REAL_STRESS_FIELD ) ;
		double bulk = (stress[0] + stress[1])/2. ;
		double shear = std::sqrt( (stress[0]-stress[1])*(stress[0]-stress[1]) + stress[0]*stress[0] + stress[1]*stress[1] )/2. ;
		out << f.getCurrentTime() << "\t" << bulk/1e6 << "\t" << shear/1e6 << "\t" << "\t" << f.getAverageField( STRAIN_FIELD, 1. )[1]*1e3 << "\t" << f.getAverageField( TENSOR_DAMAGE_FIELD, 1. )[0]*1e2 << "\t" << f.getAverageField( TENSOR_DAMAGE_FIELD, 1. )[1]*1e2 << std::endl ;
	}

	return 0 ;
}

