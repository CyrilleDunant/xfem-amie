// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../physics/materials/paste_behaviour.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../features/sample.h"
#include "../../utilities/parser.h"

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
	CommandLineParser parser("Test the logarithmic-creep material behaviour on two elements") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;


        Sample rect(nullptr, 0.04,0.04,0,0) ;
	LogarithmicCreepWithExternalParameters * behaviour = new LogarithmicCreepWithExternalParameters("young_modulus = 10e9, poisson_ratio = 0.2, imposed_deformation = 0" ) ; 
	behaviour->addMaterialLaw( new SpaceTimeDependentExternalMaterialLaw("imposed_deformation", "0.0001 t *")) ;

	rect.setBehaviour(behaviour) ;

	FeatureTree f(&rect) ;
	f.setSamplingNumber(0) ;
	int deltaTime = 1 ;
        f.setDeltaTime(1) ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
	
	std::ofstream out ;
	if(renew)
		out.open("../examples/test/test_logcreep_imposed_deformation_base", std::ios::out) ;
	else
		out.open("../examples/test/test_logcreep_imposed_deformation_current", std::ios::out) ;

	for(size_t i = 0 ; i < 5 ; i++)
	{
                f.setDeltaTime(deltaTime++) ;
		f.step() ;
		out << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, -1,1 )[1]/1e6 << "\t" << f.getAverageField( STRAIN_FIELD, -1,1 )[1]*1e3 << std::endl ;
	}

	return 0 ;
}

