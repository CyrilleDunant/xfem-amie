// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
#include "../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/materials/paste_behaviour.h"
#include "../utilities/writer/triangle_writer.h"
#include "../utilities/parser.h"
#include "../features/sample.h"

#include <fstream>
#include <omp.h>
#include <cmath>
#include <dirent.h>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

int getDelta( std::string base, std::string current, double tol )
{
	std::fstream bstream ;
	bstream.open( base.c_str(), std::ios::in ) ;
	if(!bstream.good())
	{
		std::cout << "!!! file not found: " << base << std::endl ;
		return 1 ;
	}
	std::fstream cstream ;
	cstream.open( current.c_str(), std::ios::in ) ;
	if(!cstream.good())
	{
		std::cout << "!!! file not found: " << current << std::endl ;
		return 1 ;
	}
	double bvalue = 0 ;
	double cvalue = 0 ;
	int delta = 0 ;
	while(!bstream.eof())
	{
		bstream >> bvalue ;
		if(!cstream.eof())
		{
			cstream >> cvalue ;
			if(bvalue > POINT_TOLERANCE)
			{
				if(std::abs(1.-cvalue/bvalue) > tol)
					delta++ ;
			}
			else if(std::abs(cvalue-bvalue) > tol)
				delta++ ;
		}
		else
		{
			delta++ ;
		}
	}
	while(!cstream.eof())
	{
		cstream >> cvalue ;
		delta++ ;
	}
	bstream.close() ;
	cstream.close() ;
	return delta ;
}

int main(int argc, char *argv[])
{
	CommandLineParser parser("Run all tests found in AMIE test base") ;
	parser.addFlag("--zero", false, "set tolerance to 0") ;
	parser.addValue("--tolerance", 0.01, "set tolerance for evaluation of test success" ) ;
	parser.parseCommandLine(argc, argv) ;
	double tol = std::abs(parser.getValue("--tolerance")) ;
	if(parser.getFlag("--zero"))
		tol = 0 ; 

	std::string path("../examples/test/") ;
	std::vector<std::string> exec ;
	DIR * dp ;
	struct dirent *dirp ;
	if((dp = opendir(path.c_str())) == NULL)
		std::cout << "test directory not found!" << std::endl ;

	while((dirp = readdir(dp)) != NULL)
	{
		std::string test = dirp->d_name ;
		if(test.find(".cpp") == test.size()-4 )
		{
			test.erase(test.begin(), test.begin()+5) ;
			test.erase(test.end()-4, test.end()) ;
			exec.push_back(test) ;
		}
	}

	std::cout << "cleaning existing results..." << std::endl ;
	std::system("rm ../examples/test/*_current") ;

	for(size_t i = 0 ; i < exec.size() ; i++)
	{
		timeval time0, time1 ;
		gettimeofday ( &time0, nullptr );
		std::cout << "--------------" << std::endl ;
		std::cout << "starting test " << exec[i] << std::endl ;
		std::string command = "./" + exec[i]+" 1>"+exec[i]+".out 2>"+exec[i]+".err" ;
		std::system(command.c_str()) ;
		gettimeofday ( &time1, nullptr );
		double dt = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
		std::cout << "run time: " << dt/1000000 << " seconds" << std::endl ;

		int delta = getDelta( "../examples/test/"+exec[i]+"_base", "../examples/test/"+exec[i]+"_current", tol) ;
		if(delta == 0)
			std::cout << "SUCCESS" << std::endl ;
		else
			std::cout << "FAIL: " << delta << " error(s) found" << std::endl ;
		
		
	}

	return 0 ;
}

