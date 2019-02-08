// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/stiffness.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../features/sample.h"
#include "../../features/polygonSample.h"
#include "../../features/inclusion.h"
#include "../../utilities/parser/command_line_parser.h"
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
	CommandLineParser parser("Builds a sample with a polycrystal microstructure") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addString("--output-directory","../examples/mesh/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::string dir = parser.getString("--output-directory") ;

        RectangularFeature rect(nullptr, 0.1,0.1,0,0) ;
	rect.setBehaviour(new Stiffness( 10e9, 0.2 ) ) ;
        
	std::vector<VoronoiGrain> grains ;
	grains.push_back(VoronoiGrain(new WeibullDistributedElasticStiffness( 15e9, 0.2, 0.3), 0.01, 0.6, 5, 0.55) ) ;
	grains.push_back(VoronoiGrain(new WeibullDistributedElasticStiffness( 20e9, 0.2, 0.3), 0.005, 0.3, 10, 0.75) ) ;
	grains.push_back(VoronoiGrain(new WeibullDistributedElasticStiffness( 25e9, 0.2, 0.3), 0.003, 0.1, 25) ) ;

	FeatureTree f(&rect) ;
	PSDGenerator::get2DVoronoiPolygons( &f, grains, 0, 0.003, 0.01, 32, true, 0) ;
	f.setSamplingNumber( 64 ) ;
	f.setSamplingRestriction(0.002) ;
	
	parser.setFeatureTree(&f) ;
	
	f.step() ;

	std::string name = dir+"/mesh_polycrystal_" ;
	if(renew)
		name += "base" ;
	else
		name += "current" ;

	TriangleWriter trg( name.c_str(), &f, 1. ) ;
	trg.getField( TWFT_STIFFNESS ) ;
	trg.write() ;

	return 0 ;
}

