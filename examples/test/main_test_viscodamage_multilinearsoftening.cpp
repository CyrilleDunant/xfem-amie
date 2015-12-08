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
#include "../../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
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
	CommandLineParser parser("Test a multi-linear softening damage behaviour on two elements") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::string outdir = parser.getString("--output-directory") ;


        Sample rect(nullptr, 0.01,0.01,0,0) ;
	Matrix C = PasteBehaviour(true, false, 10e9).param ;
	Point t1(0.0001, 1e6) ;
	Point t2(0.00015, 0.5e6) ;
	Point t3(0.00025, 0.1e6) ;
	std::vector<Point> tension ;
	tension.push_back(t1) ;
	tension.push_back(t2) ;
	tension.push_back(t3) ;
	std::vector<Point> compression ;
	rect.setBehaviour(new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, C, C, C*25 , new AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( tension, compression, 10e9  ), new SpaceTimeFiberBasedIsotropicLinearDamage( 0.01, 0.01, 1. ))) ;

	FeatureTree f(&rect) ;
	f.setSamplingNumber(0) ;
        f.setDeltaTime(0.001) ;
        f.setMinDeltaTime(1e-5) ;

	f.step() ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
	BoundingBoxDefinedBoundaryCondition * disp = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, 0) ;
	f.addBoundaryCondition( disp ) ;

	std::ofstream out ;
	if(renew)
		out.open(outdir+"/test_viscodamage_multilinearsoftening_base", std::ios::out) ;
	else
		out.open(outdir+"/test_viscodamage_multilinearsoftening_current", std::ios::out) ;

	for(size_t i = 0 ; i < 30 ; i++)
	{
		disp->setData( 0.0000009 + 0.0000001*i ) ;
		f.step() ;
		out << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, -1,1 )[1]/1e6 << "\t" << f.getAverageField( STRAIN_FIELD, -1,1 )[1]*1e3 << std::endl ;
	}

	return 0 ;
}

