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
#include "../../physics/damagemodels/spacetimefiberbasedfixedcrack.h"
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
	CommandLineParser parser("Test a fixed crack damage behaviour on two elements") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::string outdir = parser.getString("--output-directory") ;

	std::ofstream out ;
	if(renew)
		out.open(outdir+"/test_spacetime_fixed_crack_base", std::ios::out) ;
	else
		out.open(outdir+"/test_spacetime_fixed_crack_current", std::ios::out) ;

        RectangularFeature rect(nullptr, 0.01,0.01,0,0) ;
        SpaceTimeLimitSurfaceFractureCriterion * tension = new SpaceTimeLimitSurfaceFractureCriterion( "sqrt ( stress_1 * stres_1 + stress_2 * stress_2 )", "1e6", "", FRAME_PRINCIPAL, true ) ;
        SpaceTimeFiberBasedFixedCrack * dam = new SpaceTimeFiberBasedFixedCrack( 0.01, 1e-6, 0.99999, REAL_STRESS_FIELD, false, false ) ;

	rect.setBehaviour(new LogarithmicCreepWithExternalParameters( "young_modulus = 10e9, poisson_ratio = 0.2", tension, dam ) ) ;

	FeatureTree f(&rect) ;
	f.setSamplingNumber(0) ;
        f.setDeltaTime(0.001) ;
        f.setMinDeltaTime(1e-9) ;

	f.step() ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT_AFTER) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
	BoundingBoxDefinedBoundaryCondition * top = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, 0) ;
	BoundingBoxDefinedBoundaryCondition * right = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_XI, RIGHT_AFTER, 0) ;
	f.addBoundaryCondition( top ) ;
	top->setData( 0.0000009 ) ;

	f.step() ;

	Vector stress = f.getAverageField( REAL_STRESS_FIELD, 1. ) ;
	Vector strain = f.getAverageField( TOTAL_STRAIN_FIELD, 1. ) ;
	Vector omega = f.getAverageField( TENSOR_DAMAGE_FIELD, 1. ) ;
	out << top->getData()*1e6 << "\t" << right->getData()*1e6 << "\t" << stress[0]/1e6  << "\t" << stress[1]/1e6 << "\t" << stress[2]/1e6 << "\t" << strain[0]*1e3  << "\t" << strain[1]*1e3 << "\t" << strain[2]*1e3 << "\t" << omega[0]*1e2 << "\t" << omega[1]*1e2 << std::endl ;

	for(double i = 0. ; i < 7 ; i+=2)
	{

		top->setData( 0.000001 + 0.0000001*i ) ;

		f.step() ;

		stress = f.getAverageField( REAL_STRESS_FIELD, 1. ) ;
		strain = f.getAverageField( TOTAL_STRAIN_FIELD, 1. ) ;
		omega = f.getAverageField( TENSOR_DAMAGE_FIELD, 1. ) ;
		out << top->getData()*1e6 << "\t" << right->getData()*1e6 << "\t" << stress[0]/1e6  << "\t" << stress[1]/1e6 << "\t" << stress[2]/1e6 << "\t" << strain[0]*1e3  << "\t" << strain[1]*1e3 << "\t" << strain[2]*1e3 << "\t" << omega[0]*1e2 << "\t" << omega[1]*1e2 <<  std::endl ;
	}

//	f.removeBoundaryCondition( top ) ;
//	top->setData( 0. ) ;

	f.step() ;

	stress = f.getAverageField( REAL_STRESS_FIELD, 1. ) ;
	strain = f.getAverageField( TOTAL_STRAIN_FIELD, 1. ) ;
	omega = f.getAverageField( TENSOR_DAMAGE_FIELD, 1. ) ;
	out << top->getData()*1e6 << "\t" << right->getData()*1e6 << "\t" << stress[0]/1e6  << "\t" << stress[1]/1e6 << "\t" << stress[2]/1e6 << "\t" << strain[0]*1e3  << "\t" << strain[1]*1e3 << "\t" << strain[2]*1e3 << "\t" << omega[0]*1e2 << "\t" << omega[1]*1e2 << std::endl ;

	f.addBoundaryCondition( right ) ;

	right->setData( 0.0000009 ) ;

	f.step() ;

	stress = f.getAverageField( REAL_STRESS_FIELD, 1. ) ;
	strain = f.getAverageField( TOTAL_STRAIN_FIELD, 1. ) ;
	omega = f.getAverageField( TENSOR_DAMAGE_FIELD, 1. ) ;
	out << top->getData()*1e6 << "\t" << right->getData()*1e6 << "\t" << stress[0]/1e6  << "\t" << stress[1]/1e6 << "\t" << stress[2]/1e6 << "\t" << strain[0]*1e3  << "\t" << strain[1]*1e3 << "\t" << strain[2]*1e3 << "\t" << omega[0]*1e2 << "\t" << omega[1]*1e2 << std::endl ;

	for(double i = 0. ; i < 4 ; i+=1.5)
	{

		right->setData( 0.000001 + 0.0000001*i ) ;

		f.step() ;

		stress = f.getAverageField( REAL_STRESS_FIELD, 1. ) ;
		strain = f.getAverageField( TOTAL_STRAIN_FIELD, 1. ) ;
		omega = f.getAverageField( TENSOR_DAMAGE_FIELD, 1. ) ;
		out << top->getData()*1e6 << "\t" << right->getData()*1e6 << "\t" << stress[0]/1e6  << "\t" << stress[1]/1e6 << "\t" << stress[2]/1e6 << "\t" << strain[0]*1e3  << "\t" << strain[1]*1e3 << "\t" << strain[2]*1e3 << "\t" << omega[0]*1e2 << "\t" << omega[1]*1e2 <<  std::endl ;
	}


	return 0 ;
}

