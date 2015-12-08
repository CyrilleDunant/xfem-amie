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
#include "../utilities/font.h"
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

int getDelta( std::string base, std::string current, double tol, double ignore )
{
	std::fstream bstream ;
	bstream.open( base.c_str(), std::ios::in ) ;
	if(!bstream.good())
	{
//		std::cout << "!!! file not found: " << base << std::endl ;
		return -2 ;
	}
	std::fstream cstream ;
	cstream.open( current.c_str(), std::ios::in ) ;
	if(!cstream.good())
	{
//		std::cout << "!!! file not found: " << current << std::endl ;
		return -1 ;
	}
	double bvalue = 0 ;
	double cvalue = 0 ;
	std::string raw ;
	int delta = 0 ;
	while(!bstream.eof())
	{
		bstream >> bvalue ;
		cstream >> raw ;
		if(!bstream.eof() && !cstream.eof())
		{
                        cvalue = atof(raw.c_str()) ;
			if(std::abs(cvalue) > ignore)
			{
				if(std::abs(1.-cvalue/bvalue) > tol)
					delta++ ;
			}
			else if(std::abs(cvalue-bvalue) > ignore)
				delta++ ;
			else if( raw == "nan" || raw == "-nan")
				delta++ ;
		}
		else if(!bstream.eof())
		{
			delta++ ;
		}
	}
	while(!cstream.eof())
	{
		delta++ ;
		cstream >> raw ;
	}
	bstream.close() ;
	cstream.close() ;
	return delta ;
}

bool isDeprecated( std::string test )
{
	std::string base = "../examples/test/"+test+"_base" ;
	std::fstream in ;
	in.open(base.c_str(), std::ios::in ) ;
	std::string buffer ;
	in >> buffer ;
	return buffer == "deprecated";
}

int main(int argc, char *argv[])
{
	CommandLineParser parser("Run all tests found in AMIE test base") ;
	parser.addFlag("--zero", "set tolerance and thresholds to 0", "-O") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addValue("--tolerance", 0.01, "set relative tolerance for evaluation of test success (default: 0.01)", "-tol" ) ;
	parser.addValue("--threshold", 1e-6, "set absolute threshold below which success is not evaluated (default: 1e-6)", "-thres" ) ;
	parser.addValue("--timeout", 200, "maximum time (in seconds) spent for each test; use negative values for no time limit (default: 200s)" ) ;
	parser.addString("--match", "*", "runs only the tests matching the required string (default: runs all tests found)", "-m" ) ;
	parser.addString("--test", "*", "runs a single test with required string", "-t" ) ;
	parser.addString("--disable", "", "forces a specific test to be skipped for all users (developers only!)" ) ;
	parser.disableFeatureTreeArguments() ;
	parser.parseCommandLine(argc, argv) ;
	double tol = std::abs(parser.getValue("--tolerance")) ;
	double thr = std::abs(parser.getValue("--threshold")) ;
	double timeout = std::abs(parser.getValue("--timeout")) ;
	std::string regexp = parser.getString("--match") ;
	std::string exact = parser.getString("--test") ;
	std::string skip = parser.getString("--disable") ;
	bool renew = parser.getFlag("--renew-base") ;
	if(parser.getFlag("--zero"))
	{
		tol = 0 ; 
		thr = 0 ;
	}

	if(skip.size() > 0)
	{
		std::cout << "Warning: you are about to disable the test " << skip << std::endl ;
		std::cout << "This will affect all users. Do you wish to continue? [y/n]" << std::flush ;
		std::string buffer ;
		getline( std::cin, buffer ) ;
		if(buffer.size() == 0 || buffer == "y" || buffer == "yes" || buffer == "Y" || buffer == "YES")
		{
			std::cout << Font(BOLD, RED) << "disabling test " << skip << Font() << std::endl ;
			std::string base = "../examples/test/test_"+skip+"_base" ;
			std::fstream out ;
			out.open(base.c_str(), std::ios::out ) ;
			out << "deprecated" ;
			out.close() ;
			return 0 ;
		}
		else
		{
			std::cout << "cancelled by user. Exiting now." <<std::endl ;
			exit(0) ;
		}
	}

	if(renew)
	{
		std::cout << "Warning: you are about to renew the base of results. Do you wish to continue? [y/n]" << std::flush ;
		std::string buffer ;
		getline( std::cin, buffer ) ;
		if(buffer.size() == 0 || buffer == "y" || buffer == "yes" || buffer == "Y" || buffer == "YES")
		{
			std::cout << Font(BOLD, RED) << "overriding existing results..." << Font() << std::endl ;
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
	std::vector<std::string> files ;
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
			if(exact != std::string("*"))
			{
				std::string base = test ;
				base.erase(base.begin(), base.begin()+5) ;
				if( base == exact )
				{
					files.push_back(test) ;
					test[test.find("_")] = '/' ;
					exec.push_back(test) ;
				}
			}
			else if( regexp == "*" || test.find(regexp) != std::string::npos ) 
			{
				files.push_back(test) ;
				test[test.find("_")] = '/' ;
				exec.push_back(test) ;
			}
		}
	}

	if(!renew)
	{
		std::cout << "cleaning existing results..." << std::endl ;
		std::system("rm ../examples/test/*_current") ;
	}

	for(size_t i = 0 ; i < exec.size() ; i++)
	{
		if(isDeprecated(files[i]))
		{
			std::cout << exec[i] << Font(BOLD, BLUE) << " skipped" << Font() <<  std::endl ;
			continue ;
		}

		timeval time0, time1 ;
		gettimeofday ( &time0, nullptr );
//		std::cout << "--------------" << std::endl ;
		std::cout << exec[i] << std::flush ;
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
		std::cout << " (" << dt/1000000 << "s)" << std::flush ;

		if(!renew)
		{
			int delta = getDelta( "../examples/test/"+files[i]+"_base", "../examples/test/"+files[i]+"_current", tol, thr) ;
			if(delta == 0 && r == 0)
				std::cout << Font(BOLD, GREEN) << " SUCCESS" << Font() << std::endl ;
			else if( r == 0 && delta > 0)
				std::cout << Font(BOLD, RED) << " FAIL " << Font() << delta << " error(s) found" << std::endl ;
			else if( r == 0 && delta < 0)
				std::cout << Font(BOLD, RED) << " FAIL " << Font() << Font(BLUE) << "file not found: " << exec[i] + (r==-1 ? "_base" : "_current") << Font() <<  std::endl ;
			else
				std::cout << Font(BOLD, RED) << " FAIL " << Font() << "return value " << r << std::endl ;
		}
		else
		{
			if(r == 0)
				std::cout << Font(BOLD, GREEN) << " SUCCESS" << Font() << std::endl ;
			else
				std::cout << Font(BOLD, RED) << " FAIL " << Font() << "return value " << r << std::endl ;
		}
		
		
	}

	return 0 ;
}

