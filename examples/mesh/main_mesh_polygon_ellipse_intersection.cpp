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
	CommandLineParser parser("Builds a sample with a regular polygonal inclusion intersection an ellipsoidal inclusion") ;
	parser.addFlag("--renew-base", false, "renew the base of results") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;

        Sample rect(nullptr, 0.04,0.04,0,0) ;
	rect.setBehaviour(new Stiffness( 10e9, 0.2 ) ) ;
        
	std::valarray<Point *> vertex(7) ;
	for(size_t i = 0 ; i < vertex.size() ; i++)
		vertex[i] = new Point( 0.014*cos(0.1+2.*M_PI*i/vertex.size())-0.005, 0.012*sin(0.1+2.*M_PI*i/vertex.size())-0.004 ) ;
	PolygonalSample inc(&rect, vertex) ;
	inc.setBehaviour(new Stiffness(20e9, 0.2) ) ;

	EllipsoidalInclusion inc2(Point(0.004,0.005), Point(0.01,0.01), Point(0.007,-0.007)) ;
	inc2.setBehaviour(new Stiffness(25e9, 0.2) ) ;

	FeatureTree f(&rect) ;
	f.addFeature(&rect, &inc) ;
	f.addFeature(&inc, &inc2) ;
	f.setSamplingNumber(8) ;
	
	parser.setFeatureTree(&f) ;
	
	f.step() ;

	std::string name = "../examples/mesh/mesh_polygon_ellipse_intersection_" ;
	if(renew)
		name += "base" ;
	else
		name += "current" ;

	TriangleWriter trg( name.c_str(), &f, 1. ) ;
	trg.getField( TWFT_STIFFNESS ) ;
	trg.write() ;

	return 0 ;
}

