// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2015
//         Alain Giorla <alain.b.giorla@gmail.com>, (C) 2009-2015
//
// Copyright: See COPYING file that comes with this distribution
//
// ==== Description ====
//
// Template to create visual tests for the mesher
// Parts requiring modification are marked with an arrow "=>"
//
// ==== Recommendations ====
//
// The goal of these tests is to generate a mesh with inclusions in different conditions.
// Mechanical behaviours with imposed stresses or strains or with damage should be avoided.
// Boundary conditions should also remain as simple as possible.
// Stiffness behaviours with no external boundary conditions should be preferred for fast element generation and resolution.
// Typical run time should be a few seconds.
//
// ==== Name conventions ====
//
// Make sure the name of the TriangleWriter matches the name of the *.cpp file.
// The *.cpp file should have the following pattern: "main_mesh_*.cpp"
// The triangle file should then be "../examples/mesh/mesh_*_base" or "../examples/mesh/mesh_*_current"
// (with * the actual name of the test)
//
// ==== Integration ====
//
// 1. Move the *.cpp file in "trunk/examples/mesh"
// 2. Add a new executable in "trunk/CMakeLists.txt" along the other mesh tests
// 3. Add it to the custom target "meshTests" in "trunk/CMakeLists.txt"
// 4. Compile
// 5. Run with command line argument "--renew-base" to generate the base mesh
// 6. Run "./check_mesh" and check if the triangle files are properly generated
// 7. Add the *.cpp file and the base triangle file to the SVN repository
// 8. Commit the changes to the SVN repository (including "trunk/CMakeLists.txt")
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/stiffness.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../features/sample.h"
#include "../../utilities/parser.h"

#include <fstream>
#include <omp.h>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

int main(int argc, char *argv[])
{
	// initialize command line parser
	// => provide short description of the problem and the expected results if necessary
	CommandLineParser parser("Description of the current problem") ; 

	// command line arguments; a test should be able to run by itself without any additional command line argument
	parser.addFlag("--renew-base", false, "renew the base of results") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;

	// => define sample
        Sample rect(nullptr, 0.04,0.04,0,0) ;
	rect.setBehaviour(new Stiffness( 10e9, 0.2 ) ) ;

	// => define features here 




	// initialize FeatureTree
	// => choose default sampling number, sampling restriction, and other meshing parameters here
	FeatureTree f(&rect) ;
	f.setSamplingNumber(16) ;

	// => add features to FeatureTree here
	



	// synchronize FeatureTree with command line argument
	// this allows to set several parameters from the command line
	// notably, it can be used to verify that a geometry can be meshed consistently with several sampling numbers
	// run with "--help" to see the list of possible parameters
	parser.setFeatureTree(&f) ;
	
	// initialize elements and make first step
	f.step() ;

	// create result file
	// => modify "template" with actual name of the test
	std::string name = "../examples/mesh/mesh_template_" ;
	if(renew)
		name += "base" ;
	else
		name += "current" ;

	// initialize export
	TriangleWriter trg( name.c_str(), &f, 1. ) ;
	trg.getField( TWFT_STIFFNESS ) ;
	// => add additional fields if necessary

	// actual export
	trg.write() ;

	return 0 ;
}

