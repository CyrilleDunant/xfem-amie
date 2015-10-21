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
	CommandLineParser parser("Test the time continuity for viscoelastic materials on two elements") ;
	parser.addFlag("--renew-base", false, "renew the base of results") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;


        Sample rect(nullptr, 0.01,0.01,0,0) ;
	Matrix C = ElasticOnlyPasteBehaviour(10e9).param ;
	rect.setBehaviour(new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, C, C, C*25)) ;

	FeatureTree f(&rect) ;
	f.setSamplingNumber(0) ;
        f.setDeltaTime(0.001) ;

	f.step() ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
	BoundingBoxDefinedBoundaryCondition * disp = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, 0) ;
	f.addBoundaryCondition( disp ) ;

	std::ofstream out ;
	if(renew)
		out.open("../examples/test/test_time_continuity_base", std::ios::out) ;
	else
		out.open("../examples/test/test_time_continuity_current", std::ios::out) ;

	for(size_t i = 0 ; i < 3 ; i++)
	{
		disp->setData( 0.00000075 + 0.00000075*i ) ;
		f.step() ;
		out << f.getCurrentTime() << "\t" << f.getAverageField( STRAIN_FIELD, -1,-1 )[1]*1e3 << "\t" << f.getAverageField( STRAIN_FIELD, -1,1 )[1]*1e3 << std::endl ;
	}

	return 0 ;
}

