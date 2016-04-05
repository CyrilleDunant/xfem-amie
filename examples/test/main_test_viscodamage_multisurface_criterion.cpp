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
#include "../../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
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
	CommandLineParser parser("Test a multi-surface damage behaviour on two elements") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::string outdir = parser.getString("--output-directory") ;


        Sample rect(nullptr, 0.01,0.01,0,0) ;
	SpaceTimeMultiSurfaceFractureCriterion * crit = new SpaceTimeMultiSurfaceFractureCriterion() ;
	SpaceTimeLimitSurfaceFractureCriterion * tension = new SpaceTimeLimitSurfaceFractureCriterion( "max stress_1 stress_2", "tensile_strength", "tensile_strength", PRINCIPAL_POSITIVE ) ;
	SpaceTimeLimitSurfaceFractureCriterion * compression = new SpaceTimeLimitSurfaceFractureCriterion( "max stress_1 stress_2", "compressive_strength", "compressive_strength", PRINCIPAL_NEGATIVE ) ;
	crit->add( tension ) ;
	crit->add( compression ) ;
	SpaceTimeFiberBasedIsotropicLinearDamage * dam = new SpaceTimeFiberBasedIsotropicLinearDamage( 0.01, 1e-6, 0.99999 ) ;

	rect.setBehaviour(new LogarithmicCreepWithExternalParameters( "young_modulus = 10e9, poisson_ratio = 0.2, tensile_strength = 1e6, compressive_strength = 10e6", crit, dam ) ) ;

	FeatureTree f(&rect) ;
	f.setSamplingNumber(0) ;
        f.setDeltaTime(0.001) ;
        f.setMinDeltaTime(1e-9) ;

	f.step() ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
	BoundingBoxDefinedBoundaryCondition * disp = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, 0) ;
	f.addBoundaryCondition( disp ) ;

	std::ofstream out ;
	if(renew)
		out.open(outdir+"/test_viscodamage_multisurface_criterion_base", std::ios::out) ;
	else
		out.open(outdir+"/test_viscodamage_multisurface_criterion_current", std::ios::out) ;

	for(size_t i = 0 ; i < 10 ; i++)
	{
		disp->setData( 0.0000009 + 0.0000001*i ) ;
		f.step() ;
		out << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, 1. )[1]/1e6 << "\t" << f.getAverageField( STRAIN_FIELD, 1. )[1]*1e3 << "\t" << f.getAverageField( SCALAR_DAMAGE_FIELD, 1. )[0]*1e2 << std::endl ;
	}
	for(size_t i = 0 ; i < 20 ; i++)
	{
		disp->setData( 0.000001 - 0.000002*i ) ;
		f.step() ;
		out << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, 1. )[1]/1e6 << "\t" << f.getAverageField( STRAIN_FIELD, 1. )[1]*1e3 << "\t" << f.getAverageField( SCALAR_DAMAGE_FIELD, 1. )[0]*1e2 << std::endl ;
	}

	return 0 ;
}

