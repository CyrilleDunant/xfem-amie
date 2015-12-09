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
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addFlag("--constant","sets a constant boundary condition", "-c") ;
	parser.addFlag("--free","sets no boundary condition", "-f") ;
	parser.addFlag("--stress","sets the mechanical boundary condition in stress instead of imposed displacements", "-S") ;
	parser.addFlag("--incremental-delta-time","progressively increments the time step" ) ;
	parser.addValue("--steps", 10, "number of loading steps","-s") ;
	parser.addValue("--maximum-value", 0.001, "maximum value of the mechanical boundary condition","-m") ;
	parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	std::string test = parser.getString("--input-file") ;
	Form * behaviour = parser.getBehaviour( "--input-file" , new Stiffness(10e9, 0.2), SPACE_TWO_DIMENSIONAL ) ;
//	parser.parseConfigFile( test ) ;

	bool constant = parser.getFlag("--constant") ;
	bool free = parser.getFlag("--free") ;
	bool stress = parser.getFlag("--stress") ;
	bool renew = parser.getFlag("--renew-base") ;
	bool incr = parser.getFlag("--incremental-delta-time") ;
	size_t steps = parser.getValue("--steps") ;
	double val = parser.getValue("--maximum-value") ;
	double init = (constant ? val : 0 ) ;
	std::string outdir = parser.getString("--output-directory") ;
	test = test.substr(test.rfind('/')+1, std::string::npos) ;
	test = test.substr(0, test.length()-4) ;
	
	Sample sample(0.01,0.01,0,0) ;
	sample.setBehaviour(behaviour) ;

	FeatureTree F(&sample) ;
	
	parser.setFeatureTree(&F) ;
	F.setSamplingNumber(0) ;

	F.step() ;
	std::cout << F.getCurrentTime() << "\t" << F.getAverageField(STRAIN_FIELD, -1, 1)[1]*1e3 << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1)[1]/1e6 << std::endl ;
	double dt = F.getDeltaTime() ;

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

	std::fstream out ;
	if(renew)
	{
		out.open( outdir+"/check_behaviour_"+test+"_base", std::ios::out ) ;
	}
	else
	{
		out.open( outdir+"/check_behaviour_"+test+"_current", std::ios::out ) ;
	}

	for(size_t i = 1 ; i < steps+1 ; i++)
	{	
		if(!constant)
			up->setData(val*(double) i/(double) steps) ;
		if(incr)
			F.setDeltaTime( F.getDeltaTime()+dt ) ;
		F.step() ;
		std::cout << F.getCurrentTime()  << "\t" << F.getAverageField(STRAIN_FIELD, -1, 1)[1]*1e3 << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1)[1]/1e6 << std::endl ;
		out << F.getCurrentTime() << "\t" << F.getAverageField(STRAIN_FIELD, -1, 1)[1]*1e3 << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1)[1]/1e6 << std::endl ;
	}		

	F.getAssembly()->print() ;

	return 0 ;
}

