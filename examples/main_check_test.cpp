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

int getDelta( std::string base, std::string current, double tol, double ignore )
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
			if(std::abs(bvalue) > ignore)
			{
				if(std::abs(1.-cvalue/bvalue) > tol)
					delta++ ;
			}
			else if(std::abs(cvalue-bvalue) > ignore)
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
	parser.addFlag("--zero", false, "set tolerance and thresholds to 0") ;
	parser.addFlag("--renew-base", false, "renew the base of results") ;
	parser.addValue("--tolerance", 0.01, "set relative tolerance for evaluation of test success (default: 0.01)" ) ;
	parser.addValue("--threshold", 1e-8, "set absolute threshold below which success is not evaluated (default: 1e-8)" ) ;
	parser.addValue("--timeout", 10, "maximum time (in seconds) spent for each test; use negative values for no time limit (default: 10s)" ) ;
	parser.parseCommandLine(argc, argv) ;
	double tol = std::abs(parser.getValue("--tolerance")) ;
	double thr = std::abs(parser.getValue("--threshold")) ;
	double timeout = std::abs(parser.getValue("--timeout")) ;
	if(parser.getFlag("--zero"))
	{
		tol = 0 ; 
		thr = 0 ;
	}

	bool renew = parser.getFlag("--renew-base") ;
	if(renew)
	{
		std::cout << "Warning: you are about to renew the base of results. Do you wish to continue? [y/n]" << std::flush ;
		std::string buffer ;
		getline( std::cin, buffer ) ;
		if(buffer.size() == 0 || buffer == "y" || buffer == "yes" || buffer == "Y" || buffer == "YES")
		{
			std::cout << "overriding existing results..." <<std::endl ;
			timeout = -1 ;
		}
		else
		{
			std::cout << "cancelled by user. Exiting now." <<std::endl ;
			exit(0) ;
		}
	}


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

	if(!renew)
	{
		std::cout << "cleaning existing results..." << std::endl ;
		std::system("rm ../examples/test/*_current") ;
	}

	for(size_t i = 0 ; i < exec.size() ; i++)
	{
		timeval time0, time1 ;
		gettimeofday ( &time0, nullptr );
		std::cout << "--------------" << std::endl ;
		std::cout << "starting test " << exec[i] << std::endl ;
		std::string command ;
		if(timeout > 0)
			command = "timeout "+itoa(timeout)+" ";
		command += "./" + exec[i] ;
		for(int j = 1 ; j < argc ; j++)
			command += " " + std::string(argv[j]) ;
		command += " 1>"+exec[i]+".out 2>"+exec[i]+".err" ;
		int r = std::system(command.c_str()) ;
		gettimeofday ( &time1, nullptr );
		double dt = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
		std::cout << "run time: " << dt/1000000 << " seconds" << std::endl ;

		if(!renew)
		{
			int delta = getDelta( "../examples/test/"+exec[i]+"_base", "../examples/test/"+exec[i]+"_current", tol, thr) ;
			if(delta == 0 && r == 0)
				std::cout << "SUCCESS" << std::endl ;
			else if( r == 0)
				std::cout << "FAIL: " << delta << " error(s) found" << std::endl ;
			else
				std::cout << "FAIL (return value " << r << ")" << std::endl ;
		}
		else
		{
			if(r == 0)
				std::cout << "SUCCESS" << std::endl ;
			else
				std::cout << "FAIL (return value " << r << ")" << std::endl ;
		}
		
		
	}

	return 0 ;
}

