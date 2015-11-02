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
#include "../../utilities/parser.h"
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
	parser.addFlag("--renew-base", false, "renew the base of results") ;
	parser.addValue("--inclusions", 1000, "number of inclusions (default: 1000)") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;

        Sample rect(nullptr, 0.1,0.1,0,0) ;
	rect.setBehaviour(new Stiffness( 10e9, 0.2 ) ) ;
        
	FeatureTree f(&rect) ;
	PSDGenerator::get2DConcrete( &f, new Stiffness(20e9, 0.2), parser.getValue("--inclusions"), 0.008, 0.00001, new PSDBolomeA(), new CircularInclusionGenerator(), 100000, 0.7, nullptr, std::vector<Geometry *>(), 1) ;
	f.setSamplingNumber( 128 ) ;
	f.setSamplingRestriction(8) ;
	
	parser.setFeatureTree(&f) ;
	
	f.step() ;

	std::string name = "../examples/mesh/mesh_concrete_circle_" ;
	if(renew)
		name += "base" ;
	else
		name += "current" ;

	TriangleWriter trg( name.c_str(), &f, 1. ) ;
	trg.getField( TWFT_STIFFNESS ) ;
	trg.write() ;

	return 0 ;
}

