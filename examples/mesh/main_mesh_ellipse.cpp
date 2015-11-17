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
	CommandLineParser parser("Builds a sample with three non-intersecting ellipsoidal inclusions with aspect ratios of 0.3, 0.6 and 0.9") ;
	parser.addFlag("--renew-base", false, "renew the base of results") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;

        Sample rect(nullptr, 0.04,0.04,0,0) ;
	rect.setBehaviour(new Stiffness( 10e9, 0.2 ) ) ;
        
	EllipsoidalInclusion ell9( Point(0,0.01), Point(0.0075,0), Point(0,0.0075*0.9)) ;
	ell9.setBehaviour(new Stiffness(20e9, 0.2) ) ;

	EllipsoidalInclusion ell6( Point(-0.01,-0.005), Point(0,0.01), Point(0.01*0.6,0)) ;
	ell6.setBehaviour(new Stiffness(20e9, 0.2) ) ;

	EllipsoidalInclusion ell3( Point(0.01,-0.005), Point(0,0.01), Point(0.01*0.3,0)) ;
	ell3.setBehaviour(new Stiffness(20e9, 0.2) ) ;

	FeatureTree f(&rect) ;
	f.addFeature(&rect, &ell9) ;
	f.addFeature(&rect, &ell6) ;
	f.addFeature(&rect, &ell3) ;
	f.setSamplingNumber(8) ;
	
	parser.setFeatureTree(&f) ;
	
	f.step() ;

	std::string name = "../examples/mesh/mesh_ellipse_" ;
	if(renew)
		name += "base" ;
	else
		name += "current" ;

	TriangleWriter trg( name.c_str(), &f, 1. ) ;
	trg.getField( TWFT_STIFFNESS ) ;
	trg.write() ;

	return 0 ;
}

