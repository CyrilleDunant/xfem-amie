// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../physics/stiffness.h"
#include "../utilities/parser.h"
#include "../features/sample.h"

#include <fstream>
#include <ostream>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <dirent.h>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

int main(int argc, char *argv[])
{
	CommandLineParser parser("Makes a tensile test on a 2 elements sample at a constant imposed displacement rate") ;
	parser.addArgument("file_name", "", "path to a *.ini file containing the behaviour of the sample" ) ;
	parser.addValue("--steps", 10, "number of loading steps","-s") ;
	parser.addValue("--maximum-value", 0.001, "maximum value of the mechanical boundary condition","-m") ;
	parser.addFlag("--constant","sets a constant boundary condition", "-c") ;
	parser.addFlag("--free","sets no boundary condition", "-f") ;
	parser.addFlag("--stress","sets the mechanical boundary condition in stress instead of imposed displacements", "-S") ;
	parser.parseCommandLine(argc, argv) ;

	Form * behaviour = parser.getBehaviour( "file_name" , new Stiffness(10e9, 0.2), SPACE_TWO_DIMENSIONAL ) ;
	size_t steps = parser.getValue("--steps") ;
	double val = parser.getValue("--maximum-value") ;
	bool constant = parser.getFlag("--constant") ;
	bool free = parser.getFlag("--free") ;
	bool stress = parser.getFlag("--stress") ;
	double init = (constant ? val : 0 ) ;
	
	Sample sample(0.01,0.01,0,0) ;
	sample.setBehaviour(behaviour) ;

	FeatureTree F(&sample) ;
	
	parser.setFeatureTree(&F) ;
	F.setSamplingNumber(1) ;

	F.step() ;
	std::cout << F.getCurrentTime() << "\t" << 0 << "\t" << F.getAverageField(STRAIN_FIELD, -1, 1)[1] << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1)[1] << std::endl ;

	BoundingBoxDefinedBoundaryCondition * up ;	
	if(F.getOrder() < CONSTANT_TIME_LINEAR)
	{
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM ) ) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT ) ) ;
		if(stress)
			up = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, init) ;
		else
			up = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP, init) ;
		if(!free)
			F.addBoundaryCondition(up) ;
	}
	else
	{
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER ) ) ;
		if(stress)
			up = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, init) ;
		else
			up = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP_AFTER, init) ;
		if(!free)
			F.addBoundaryCondition(up) ;
	}

	
	for(size_t i = 1 ; i < steps+1 ; i++)
	{	
		if(!constant)
			up->setData(val*(double) i/(double) steps) ;
		F.step() ;
		std::cout << F.getCurrentTime() << "\t" << up->getData() << "\t" << F.getAverageField(STRAIN_FIELD, -1, 1)[1] << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1)[1] << std::endl ;
	}		


	return 0 ;
}

