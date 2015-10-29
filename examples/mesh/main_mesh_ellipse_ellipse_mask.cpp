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
	CommandLineParser parser("Builds a sample with one ellipsoidal inclusion masking another") ;
	parser.addFlag("--renew-base", false, "renew the base of results") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;

        Sample rect(nullptr, 0.04,0.04,0,0) ;
	rect.setBehaviour(new Stiffness( 10e9, 0.2 ) ) ;
        
	EllipsoidalInclusion inc(Point(0.00,0.00), Point(0.015,0),Point(0,0.015*0.7)) ;
	inc.setBehaviour(new Stiffness(20e9, 0.2) ) ;

	EllipsoidalInclusion inc2(Point(-0.005,-0.005), Point(0.012*0.8,0),Point(0,0.012))  ;
	inc2.setBehaviour(new Stiffness(30e9, 0.2) ) ;
	inc2.addToMask(&inc) ;

	FeatureTree f(&rect) ;
	f.addFeature(&rect, &inc) ;
	f.addFeature(&inc, &inc2) ;
	f.setSamplingNumber(16) ;
	
	parser.setFeatureTree(&f) ;
	
	f.step() ;

	std::string name = "../examples/mesh/mesh_ellipse_ellipse_mask_" ;
	if(renew)
		name += "base" ;
	else
		name += "current" ;

	TriangleWriter trg( name.c_str(), &f, 1. ) ;
	trg.getField( TWFT_STIFFNESS ) ;
	trg.write() ;

	return 0 ;
}

