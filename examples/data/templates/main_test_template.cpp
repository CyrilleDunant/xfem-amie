// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2015
//         Alain Giorla <alain.b.giorla@gmail.com>, (C) 2009-2015
//
// Copyright: See COPYING file that comes with this distribution
//
// ==== Description ====
//
// Template to create unit tests
// Parts requiring modification are marked with an arrow "=>"
//
// ==== Recommendations ====
//
// The goal of these tests is to generate a small test to verify that a feature from AMIE works as intended.
// Sample geometry, mechanical behaviours and boundary conditions should be as simple as possible.
// Typical run time should be a few seconds.
//
// ==== Name conventions ====
//
// Make sure the name of the output matches the name of the *.cpp file.
// The *.cpp file should have the following pattern: "main_test_*.cpp"
// The triangle file should then be "../examples/test/test_*_base" or "../examples/test/test_*_current"
// (with * the actual name of the test)
//
// ==== Integration ====
//
// 1. Move the *.cpp file in "trunk/examples/test"
// 2. Add a new executable in "trunk/CMakeLists.txt" along the other unit tests
// 3. Add it to the custom target "unitTests" in "trunk/CMakeLists.txt"
// 4. Compile
// 5. Run with command line argument "--renew-base" to generate the base results
// 6. Run "./check_test" and check if the output files are properly generated
// 7. Add the *.cpp file and the base output file to the SVN repository
// 8. Commit the changes to the SVN repository (including "trunk/CMakeLists.txt")
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/stiffness.h"
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
	// => provide short description of the problem
	CommandLineParser parser("Description of the current problem") ; 

	// command line arguments; a test should be able to run by itself without any additional command line argument
	// furthermore, the test should always be run on the same mesh to provide consistent results
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.disableFeatureTreeArguments() ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;

	// => define sample
        Sample rect(nullptr, 0.04,0.04,0,0) ;
	rect.setBehaviour(new Stiffness( 10e9, 0.2 ) ) ;

	// => define features here 




	// initialize FeatureTree
	// => choose default sampling number, sampling restriction, and other meshing parameters here
	FeatureTree f(&rect) ;

	// => add features, boundary conditions, etc to FeatureTree here
	



	// create result file
	// => modify "template" with actual name of the test
	std::string name = "../examples/test/test_template_" ;
	if(renew)
		name += "base" ;
	else
		name += "current" ;

	// initialize export
	std::ofstream out ;
	out.open( name.c_str(), std::ios::out ) ;

	// => make as many steps as necessary
	for(size_t i = 0 ; i < 1 ; i++)
	{
		f.step() ;

		// => export necessary data here
		out << f.getCurrentTime() << std::endl ;
	}
	
	return 0 ;
}

