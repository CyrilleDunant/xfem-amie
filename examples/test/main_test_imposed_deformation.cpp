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
	CommandLineParser parser("Test the stiffness with imposed deformation material behaviour on two elements") ;
	parser.addFlag("--renew-base", false, "renew the base of results") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;


        Sample rect(nullptr, 0.04,0.04,0,0) ;
	rect.setBehaviour(new StiffnessWithImposedDeformation( 10e9, 0.2, 0.001 ) ) ;

	FeatureTree f(&rect) ;
	f.setSamplingNumber(0) ;
	
	f.step() ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM ) ) ;

	std::ofstream out ;
	if(renew)
		out.open("../examples/test/test_imposed_deformation_base", std::ios::out) ;
	else
		out.open("../examples/test/test_imposed_deformation_current", std::ios::out) ;

	f.step() ;
	out << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, -1,1 )[1] << "\t" << f.getAverageField( STRAIN_FIELD, -1,1 )[1] << std::endl ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP, -1e6 ) ) ;
	f.step() ;
	out << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, -1,1 )[1] << "\t" << f.getAverageField( STRAIN_FIELD, -1,1 )[1] << std::endl ;

	return 0 ;
}

