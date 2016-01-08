// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../features/sample.h"
#include "../../features/polygonSample.h"
#include "../../physics/stiffness.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../physics/material_laws/mechanical_material_laws.h"
#include "../../physics/material_laws/material_laws.h"
#include "../../utilities/parser.h"
#include "../../utilities/itoa.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../geometry/level_set.h" 
#include "../../utilities/mineral.h" 


#include <dirent.h>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#include <sys/time.h>


using namespace Amie ;


int main( int argc, char *argv[] )
{
	timeval time0, time1 ;
	gettimeofday ( &time0, nullptr );

	CommandLineParser parser("Runs a series of elastic test on a rock sample defined in a *.sci file", false, false) ;
	parser.addArgument("rock","granite_1","name of the rock to be tested") ;
	parser.addString("--directory","/home/ag3/Code/denisov/","path to the mineral and rock database", "-D") ;
	parser.addValue("--seed", 20, "random seed for the microstructure generation", "-s") ;
	parser.addValue("--show-microstructure", -1, "index of the microstructure to print", "-m") ;
	parser.parseCommandLine( argc, argv ) ;
	
	return 0 ;
}
