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
	CommandLineParser parser("Builds a sample with two polygonal inclusions sharing a common edge") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addString("--output-directory","../examples/mesh/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::string dir = parser.getString("--output-directory") ;

        RectangularFeature rect(nullptr, 0.04,0.04,0,0) ;
	rect.setBehaviour(new Stiffness( 10e9, 0.2 ) ) ;
        
	std::valarray<Point *> vertex(4) ;
	vertex[0] = new Point(0,-0.015) ;
	vertex[1] = new Point(0,0.015) ;
	vertex[2] = new Point(-0.015,0.01) ;
	vertex[3] = new Point(-0.01,-0.012) ;
	PolygonalSample inc(&rect, vertex) ;
	inc.setBehaviour(new Stiffness(20e9, 0.2) ) ;

	std::valarray<Point *> vertex2(4) ;
	vertex2[0] = new Point(0,-0.015) ;
	vertex2[1] = new Point(0,0.015) ;
	vertex2[2] = new Point(0.01,0.012) ;
	vertex2[3] = new Point(0.015,-0.01) ;
	PolygonalSample inc2(nullptr, vertex2) ;
	inc2.setBehaviour(new Stiffness(25e9, 0.2) ) ;

	FeatureTree f(&rect) ;
	f.addFeature(&rect, &inc) ;
	f.addFeature(&rect, &inc2) ;
	f.setSamplingNumber(8) ;
	
	parser.setFeatureTree(&f) ;
	
	f.step() ;

	std::string name = dir+"/mesh_polygon_polygon_edge_" ;
	if(renew)
		name += "base" ;
	else
		name += "current" ;

	TriangleWriter trg( name.c_str(), &f, 1. ) ;
	trg.getField( TWFT_STIFFNESS ) ;
	trg.write() ;

	return 0 ;
}

