// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/stiffness.h"
#include "../../physics/viscoelasticity_and_fracture.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../utilities/writer/triangle_writer.h"
//#include "../../physics/damagemodels/spacetimefiberbasedtensileshearlineardamage.h"
#include "../../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../../physics/damagemodels/spacetimefiberbasedplasticstrain.h"
//#include "../../physics/fracturecriteria/spacetimemultisurfacefracturecriterion.h"
#include "../../physics/fracturecriteria/spacetimelimitsurfacefracturecriterion.h"
#include "../../physics/fracturecriteria/maxstrain.h"
#include "../../features/sample.h"
#include "../../features/polygonSample.h"
#include "../../features/inclusion.h"
#include "../../utilities/parser/parser.h"
#include "../../utilities/parser/command_line_parser.h"
#include "../../utilities/postprocessor.h"
#include "../../utilities/itoa.h"

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
	CommandLineParser parser("Makes a tensile test on a 2 elements sample at a constant imposed displacement rate") ;
	parser.parseCommandLine(argc, argv) ;
	
	Sample sample(0.01,0.01,0,0) ;
	FeatureTree F(&sample) ;
	F.setSamplingNumber(0) ;

	Matrix C = Tensor::cauchyGreen( 10e9, 0.2, SPACE_TWO_DIMENSIONAL, PLANE_STRAIN, YOUNG_POISSON ) ;
	DamageModel * old = new SpaceTimeFiberBasedIsotropicLinearDamage( 0.01, 1e-9, 1. ) ;
	DamageModel * dam = new SpaceTimeFiberBasedPlasticStrain( 0.00001, 1e-9, 0.001 ) ;
	FractureCriterion * crit = new SpaceTimeLimitSurfaceFractureCriterion( "max ( stress_1 ) ( stress_2 )", "2e6", "", ALL  ) ; 
	ViscoelasticityAndFracture * frac = new ViscoelasticityAndFracture( GENERALIZED_KELVIN_VOIGT, C, C*2, C*10, crit, dam ) ;
	LogarithmicCreepWithExternalParameters * creep = new LogarithmicCreepWithExternalParameters("young_modulus = 10e9, poisson_ratio = 0.2, creep_modulus = 20e9, creep_characteristic_time = 1, tensile_strain = 0.0001, fracture_energy = 1500, plastic_orientation = 0.23", crit, dam ) ;
	Viscoelasticity * visc = new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, C, C*2, C*10 ) ;
	sample.setBehaviour( frac ) ;

	F.setDeltaTime(1.) ;
	F.setMaxIterationsPerStep( 10000 ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER ) ) ;
	BoundingBoxDefinedBoundaryCondition * up ;	
	up = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP_AFTER, 0) ;
	F.addBoundaryCondition( up ) ;

	F.step() ;

	for(size_t i = 0 ; i < 100 ; i++)
	{
		up->setData( 0.0000001*(i+1)) ;
		F.step() ;
		std::cout << F.getCurrentTime()  << "\t" << F.getAverageField(STRAIN_FIELD)[1] << "\t" << F.getAverageField(REAL_STRESS_FIELD)[1]/1e6 << std::endl ;		
	}


	return 0 ;
}

