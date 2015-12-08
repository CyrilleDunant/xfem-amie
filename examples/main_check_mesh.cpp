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
#include <regex>
#include <dirent.h>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

int main(int argc, char *argv[])
{
	CommandLineParser parser("Run all tests found in the mesh test database") ;
	parser.disableFeatureTreeArguments() ;
	parser.addFlag("--viewer", "open viewer to compare the meshes for each test", "-V") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addValue("--timeout", 10, "maximum time (in seconds) spent for each test; use negative values for no time limit (default: 10s)" ) ;
	parser.addString("--viewer-path", "viewer", "path to AMIE viewer" ) ;
	parser.addString("--match", "*", "runs only the tests matching the required string (default: runs all tests found)", "-m" ) ;
	parser.addString("--test", "*", "runs a single test with required string" , "-t") ;
	parser.addString("--amie-build-directory","./","directory of the build of the main AMIE distribution", "-A") ;
	parser.addString("--output-directory","../examples/mesh/","directory where the results are stored", "-D") ;
	parser.addString("--base-directory","../examples/mesh/","directory where the results are stored", "-B") ;
	parser.addValue("--set-sampling-number", 16, "set the number of points on the edge of the sample" ) ;
	parser.addValue("--set-sampling-restriction", 0, "set the number of mesh points below which small inclusions are not meshed" ) ;
	parser.parseCommandLine(argc, argv) ;
	double timeout = std::abs(parser.getValue("--timeout")) ;
	bool renew = parser.getFlag("--renew-base") ;
	bool compare = parser.getFlag("--viewer") ;
	std::string viewer = parser.getString("--viewer-path") ;
	std::string regexp = parser.getString("--match") ;
	std::string exact = parser.getString("--test") ;
	std::string dir = parser.getString("--amie-build-directory") ;
	std::string outdir = parser.getString("--output-directory") ;
	std::string basedir = parser.getString("--base-directory") ;
	std::string meshdir = dir+"../examples/mesh/" ;

	if(renew)
	{
		std::cout << "Warning: you are about to renew the base of results. Do you wish to continue? [y/n]" << std::flush ;
		std::string buffer ;
		getline( std::cin, buffer ) ;
		if(buffer.size() == 0 || buffer == "y" || buffer == "yes" || buffer == "Y" || buffer == "YES")
		{
			std::cout << Font( BOLD, RED ) << "overriding existing results..." << Font() << std::endl ;
			timeout = -1 ;
		}
		else
		{
			std::cout << "cancelled by user. Exiting now." <<std::endl ;
			exit(0) ;
		}
	}


	std::string path(meshdir) ;
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
					files.push_back("/"+test) ;
					test[test.find("_")] = '/' ;
					exec.push_back(test) ;
				}
			}
			else if( regexp == "*" || test.find(regexp) != std::string::npos ) 
			{
				files.push_back("/"+test) ;
				test[test.find("_")] = '/' ;
				exec.push_back(test) ;
			}
		}
	}

	if(!renew)
	{
		std::cout << "cleaning existing results..." << std::endl ;
		std::string rm = "rm "+outdir+"/*_current" ;
		std::system(rm.c_str()) ;
	}

	for(size_t i = 0 ; i < exec.size() ; i++)
	{
		timeval time0, time1 ;
		gettimeofday ( &time0, nullptr );
//		std::cout << "--------------" << std::endl ;
		std::cout << exec[i] << std::flush ;
		std::string command ;
		if(timeout > 0)
			command = "timeout "+itoa(timeout)+" ";
		command += dir + "/" + exec[i] ;
		for(int j = 1 ; j < argc ; j++)
			command += " " + std::string(argv[j]) ;
		command += " 1>"+dir+"/"+exec[i]+".out 2>"+dir+"/"+exec[i]+".err" ;
		int r = std::system(command.c_str()) ;
		gettimeofday ( &time1, nullptr );
		double dt = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
		std::cout << " (" << dt/1000000 << "s)" << std::flush ;

		if(!renew)
		{
			if(compare)
			{
				std::string viewerCommand = viewer + " "+basedir+files[i]+"_base | "+viewer+" "+outdir+files[i]+"_current" ;
				std::cout << viewerCommand << std::endl ;
				std::system( viewerCommand.c_str()  ) ;
			}

			if(r == 0)
				std::cout << Font(BOLD, GREEN) << " DONE" << Font() << std::endl ;
			else
				std::cout << Font(BOLD, RED) << " FAIL " << Font() << "return value " << r << std::endl ;
		}
		else
		{

			if(r == 0)
				std::cout << Font(BOLD, GREEN) << " DONE" << Font() << std::endl ;
			else
				std::cout << Font(BOLD, RED) << " FAIL " << Font() << "return value " << r << std::endl ;
		}
		
		
	}

	return 0 ;
}

