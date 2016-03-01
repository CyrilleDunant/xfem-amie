// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/stiffness.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../features/sample.h"
#include "../../features/polygonSample.h"
#include "../../features/inclusion.h"
#include "../../utilities/parser.h"
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
	CommandLineParser parser("Builds a sample with a concrete-like microstructure with circular aggregates") ;
	parser.addValue("--inclusions", 100, "number of inclusions (default: 100)") ;
	parser.parseCommandLine(argc, argv) ;

        Sample rect(nullptr, 0.1,0.1,0,0) ;
//	rect.setBehaviour(new LogarithmicCreepWithExternalParameters( "young_modulus = 10e9, poisson_ratio = 0.2, density = 2.2" ) ) ;
	rect.setBehaviour(new Stiffness( 10e9, 0.2 ) ) ;
        
	FeatureTree f(&rect) ;
//        Form * agg = new LogarithmicCreepWithExternalParameters( "young_modulus = 20e9, poisson_ratio = 0.2, density = 2.7") ;
        Form * agg = new Stiffness( 20e9, 0.2 ) ;
	std::vector<Feature *> aggs = PSDGenerator::get2DConcrete( &f, agg , parser.getValue("--inclusions"), 0.008, 0.00001, new PSDBolomeA(), new CircularInclusionGenerator(), 100000, 0.7, nullptr, std::vector<Geometry *>(), 1) ;
	f.setSamplingNumber( 24 ) ;
	f.setSamplingRestriction(0.002) ;
	f.setDeltaTime(1) ;
        f.setOrder( LINEAR ) ;	

	parser.setFeatureTree(&f) ;
	
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER ) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, 0.1 ) ) ;

	f.step() ;

	std::vector<Geometry *> geoms ;
        for(size_t i = 0 ; i < aggs.size() ; i++) { geoms.push_back( dynamic_cast<Geometry *>( aggs[i] ) ) ; }
        unsigned int aggIndex = f.get2DMesh()->generateCache( geoms ) ;
        std::vector<unsigned int> aggIndexes ; aggIndexes.push_back( aggIndex ) ;
        unsigned int cemIndex = f.get2DMesh()->generateCacheOut( aggIndexes ) ;

        MacroscopicStrainPostProcessor macro(0,0,1) ;
	Vector maxStr = macro.postProcess( &f ) ;         

        std::cout << maxStr[0] << " " << maxStr[1] << " " << maxStr[2] << std::endl ;

	return 0 ;
}

