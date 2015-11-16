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
	CommandLineParser parser("Test the viscoelasticity with imposed deformation material behaviour on two elements") ;
	parser.addFlag("--renew-base", false, "renew the base of results") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;


        Sample rect(nullptr, 0.04,0.04,0,0) ;
	Matrix M = Tensor::cauchyGreen(10e9, 0.2, true, SPACE_TWO_DIMENSIONAL, PLANE_STRESS) ;
	Vector v(3) ; v[0]=0.001 ; v[1]=0.001 ;
	rect.setBehaviour(new ViscoelasticityAndImposedDeformation( GENERALIZED_KELVIN_VOIGT, M, M, M*5, v ) ) ;
//	rect.setBehaviour(new LogarithmicCreepWithExternalParameters( "young_modulus = 10e9, poisson_ratio = 0.2, imposed_deformation = 0.001" ) ) ;

	FeatureTree f(&rect) ;
	f.setOrder(LINEAR_TIME_LINEAR) ;
	f.setSamplingNumber(0) ;
        f.setDeltaTime(1) ;
	
	f.step() ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;

	std::ofstream out ;
	if(renew)
		out.open("../examples/test/test_viscoelasticity_imposed_deformation_base", std::ios::out) ;
	else
		out.open("../examples/test/test_viscoelasticity_imposed_deformation_current", std::ios::out) ;

	f.step() ;
	out << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, -1,1 )[1]/1e6 << "\t" << f.getAverageField( STRAIN_FIELD, -1,1 )[1]*1e3 << std::endl ;

	f.setDeltaTime(2) ;
	f.step() ;
	out << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, -1,1 )[1]/1e6 << "\t" << f.getAverageField( STRAIN_FIELD, -1,1 )[1]*1e3 << std::endl ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, -1e6 ) ) ;
	f.step() ;
	out << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, -1,1 )[1]/1e6 << "\t" << f.getAverageField( STRAIN_FIELD, -1,1 )[1]*1e3 << std::endl ;

	f.setDeltaTime(3) ;
	f.step() ;
	out << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, -1,1 )[1]/1e6 << "\t" << f.getAverageField( STRAIN_FIELD, -1,1 )[1]*1e3 << std::endl ;

	return 0 ;
}

