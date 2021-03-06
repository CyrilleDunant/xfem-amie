// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/configuration.h"
#include "../utilities/parser/config_parser.h"
#include "../utilities/parser/command_line_parser.h"
#include "../utilities/postprocessor.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#include <dirent.h>

using namespace Amie ;

std::vector<std::string> getFiles( std::string path, std::string end)
{
    std::vector<std::string> ret ;

    DIR * dp ;
    struct dirent *dirp ;
    if((dp = opendir(path.c_str())) == NULL)
    {
        std::cout << "could not open directory " << path << std::endl ; 
        exit(0) ;
    }

    while((dirp = readdir(dp)) != NULL)
    {
        std::string test = dirp->d_name ;
        if(test.find(end) == test.size()-end.size() && test.size() > end.size())
            ret.push_back( path+"/"+test ) ;
    }

    return ret ;
}

int main(int argc, char *argv[])
{

    timeval time0, time1 ;
    gettimeofday ( &time0, nullptr );

    CommandLineParser parser("Run a 2D AMIE simulation using the configuration parameters found in a *.ini file", true, true) ;
    parser.addFlag( "--print-microstructure", "print only the mesh and values of mechanical properties", "-M") ;
    parser.addFlag( "--print-configuration-tree", "print only the configuration tree after parsing", "-ct") ;
    parser.addFlag( "--browse", "call the explorer to find the *.ini file", "-b") ;
    parser.addString( "--directory", std::string() , "directory to read additional input and write output", "-D") ;
    parser.addArgument( "file_name", "*.ini", "relative path to *.ini file to run") ;
    parser.parseCommandLine( argc, argv ) ;
    bool explorer = parser.getFlag("--browse") ;
    bool printMicrostructureOnly = parser.getFlag("--print-microstructure") ;
    bool printConfigTree = parser.getFlag("--print-configuration-tree") ;
    std::string file = parser.getStringArgument(0) ;
    std::string path = parser.getString("--directory") ;

    if(file == "*.ini" && !explorer)
    {
        std::string p = "." ;
        if(path.size() > 0)
            p = path ;
        std::vector<std::string> potentials = getFiles(p,".ini");
        if(potentials.size() == 0)
        {
           std::cout << "no *.ini file defined, nothing to do..." << std::endl ;
           exit(0) ;
        }
        if(potentials.size() == 1)
        {
           std::cout << "using file " << potentials[0] << " automatically detected in working directory" << std::endl ;
           file = potentials[0] ;
        }
        else
        {
            std::cout << "found multiple files in working directory:" << std::endl ;
            for(size_t i = 0 ; i < potentials.size() ; i++)
                std::cout << i << " " << potentials[i] << std::endl ;
            std::cout << "choose *.ini file (0-" << potentials.size()-1 << "; \"x\" to open explorer window)?" << std::flush ;
            std::string buffer ;
            getline( std::cin, buffer ) ;
            if(buffer == "x")
            {
                explorer = true ;
            }
            else
            {
                bool found = false ;
                for(size_t i = 0 ; i < potentials.size() ; i++)
                {
                    if(path+"/"+buffer == potentials[i])
                    {
                        found = true ;
                        file = potentials[i] ;
                    }
                }
                if(!found)
                {
                    size_t index = atof( buffer.c_str() ) ;
                    if(index >= potentials.size() || !FunctionParserHelper::isNumeral( buffer[0] ))
                    {
                        std::cout << "no *.ini file defined, nothing to do..." << std::endl ;
                        exit(0) ;
                    } 
                    file = potentials[index] ;
                }
                std::cout << "using file " << file << std::endl ;
            }
        }
    }

    if(explorer)
    {
        std::string command = "../input_manager/input_manager" ;
        for( int i = 1 ; i < argc ; i++)
        {
            if(std::string(argv[i]) != "--browse")
                command += " "+std::string(argv[i]) ;
        }
        std::system( command.c_str() ) ;
        return 0 ;
    }


    ConfigTreeItem * define = parser.getLocalConfiguration() ;
    std::map<std::string, std::string> direct = parser.getDirectConfiguration() ;
    std::vector<std::string> flags = parser.getActiveFlags() ;

    ConfigTreeItem * problem = ConfigParser::readFile(file, define, true, true, flags, path ) ;
    problem->configure( direct ) ;

    if(printConfigTree)
    {
        problem->printTree() ;
        exit(0) ;
    }

    FeatureTree F(problem->getChild("sample")->getSample( )) ;
    if(problem->hasChildFromFullLabel("sample.sampling_factor"))
        F.setSamplingFactor( F.getFeature(0), problem->getData("sample.sampling_factor", 1.) ) ;
    if(problem->hasChildFromFullLabel("sample.sampler"))
        F.setSampler( F.getFeature(0), problem->getChildFromFullLabel("sample.sampler")->getSampler() ) ;
    F.setDiscretizationParameters(problem->getChild("discretization")) ;
    Vector instants = F.setSteppingParameters(problem->getChild("stepping")) ;
    InclusionFamilyTree * allFeatures = nullptr ;
    if(problem->hasChild("inclusions"))
        allFeatures = problem->getChild("inclusions")->makeInclusionFamilyTree( &F ) ;
    else
        allFeatures = new InclusionFamilyTree( ) ;
    allFeatures->concatenate() ;

    if(problem->hasChild("cracks"))
    {
        if(F.getOrder() >= LINEAR_TIME_LINEAR)
        {
            std::cout << "cannot use XFEM cracks with space-time finite elements, exiting now!" << std::endl;
            return 0;
        }
        std::vector<BranchedCrack *> cracks = problem->getChild("cracks")->getAllCracks();
        for(size_t i = 0; i < cracks.size(); i++)
            F.addFeature( F.getFeature(0), cracks[i] );
    }

    parser.setFeatureTree( &F ) ;

    std::vector<BoundaryCondition * > bc_prestep;
    if(problem->hasChild("prestep"))
        bc_prestep = problem->getChild("prestep")->getAllBoundaryConditions(&F) ;

    F.step() ;

    for (size_t i = 0; i < bc_prestep.size(); i++)
        F.removeBoundaryCondition(bc_prestep[i]);

    std::vector<unsigned int> cacheIndex ;
    std::vector<unsigned int> aggCacheIndex ;
    std::cout << "generating cache for inclusion family 0/" << allFeatures->features.size() ;
    cacheIndex.push_back( 0 ) ;
    for(size_t i = 0 ; i < allFeatures->features.size() ; i++)
    {
        std::cout << "\rgenerating cache for inclusion family " << i+1 << "/" << allFeatures->features.size() ;
        cacheIndex.push_back( F.get2DMesh()->generateCache( allFeatures->getFeaturesAsGeometry(i) ) ) ;
        aggCacheIndex.push_back( cacheIndex[cacheIndex.size()-1] ) ;
    }
    std::cout << "generating cache for inclusion family 0/" << allFeatures->features.size() ;
    cacheIndex[0] = F.get2DMesh()->generateCacheOut( aggCacheIndex ) ;
    std::cout << "... done" << std::endl ;
    for(size_t i = 0 ; i < cacheIndex.size() ; i++)
        std::cout << "inclusion family " << i << " covering surface " << F.get2DMesh()->getArea( cacheIndex[i] ) << std::endl ;

    if(printMicrostructureOnly)
    {
        TriangleWriter trg(problem->getStringData("export.file_name","2d_composite_microstructure"), &F, 1 ) ;
        trg.getField( TWFT_STIFFNESS ) ;
        trg.write() ;
        exit(0) ;
    }

    MultiTriangleWriter * trg = nullptr ;
    if(problem->hasChildFromFullLabel("export.file_name"))
    {
        std::string trgFileName = problem->getStringData("export.file_name","file_not_found") ;
        std::string headerFileName = trgFileName + "_header" ;
        trg = new MultiTriangleWriter( headerFileName, trgFileName, &F, 1.) ;
    }
    std::vector<PostProcessor *> posts ; std::string postProcessFile = "file_not_found" ;
    if(problem->hasChild("post_processor"))
    {
        posts = problem->getChild("post_processor")->getAllPostProcessors(cacheIndex) ;
        postProcessFile = problem->getStringData("post_processor.file_name","file_not_found") ;
    }

        if(problem->hasChild("output"))
            problem->getChild("output")->writeOutput(&F, 0, instants.size(), cacheIndex, flags) ;
        if(problem->hasChild("export"))
            problem->getChild("export")->exportSvgTriangles(trg, &F, 0, instants.size(), flags) ;

    if(problem->hasChild("boundary_conditions"))
        std::vector<BoundaryCondition * > bc = problem->getChild("boundary_conditions")->getAllBoundaryConditions(&F) ;

    bool next = true ;
    for(size_t i = 1 ; i < instants.size() ; i++)
    {
        if(next)
            F.setDeltaTime( instants[i]-instants[i-1] ) ;

        next = F.step() ;
//		F.getAssembly()->print() ;
//		exit(0) ;
        if(problem->hasChild("output"))
            problem->getChild("output")->writeOutput(&F, i, instants.size(), cacheIndex, flags) ;
        if(problem->hasChild("export"))
            problem->getChild("export")->exportSvgTriangles(trg, &F, i, instants.size(), flags) ;
        if(problem->hasChild("post_processor") && postProcessFile != std::string("file_not_found"))
            PostProcessor::write( postProcessFile, &F, posts, true, i==1 ) ;
    }

    gettimeofday ( &time1, nullptr );
    double delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
    std::cout << "problem solved in " << delta/1000000 << " seconds" << std::endl ;

    parser.sendEmail("AMIE - Execution finished", parser.getCommandLine(), problem->getStringData("output.file_name",std::string()) ) ;

    return 0 ;
}
