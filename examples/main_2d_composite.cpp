// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/configuration.h"
#include "../utilities/parser.h"

#ifdef HAVE_OMP
#include <omp>
#endif
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>

using namespace Amie ;

int main(int argc, char *argv[])
{

    timeval time0, time1 ;
    gettimeofday ( &time0, nullptr );

    CommandLineParser parser("Run a 2D AMIE simulation using the configuration parameters found in a *.ini file", true, true) ;
    parser.addFlag( "--print-microstructure", false , "print only the mesh and values of mechanical properties") ;
    parser.addFlag( "--print-configuration-tree", false , "print only the configuration tree after parsing") ;
    parser.addArgument( "file_name", "../examples/data/composite/test_2d_composite.ini", "relative path to *.ini file to run") ;
    parser.parseCommandLine( argc, argv ) ;
    bool printMicrostructureOnly = parser.getFlag("--print-microstructure") ;
    bool printConfigTree = parser.getFlag("--print-configuration-tree") ;
    std::string file = parser.getStringArgument(0) ;

    ConfigTreeItem * define = parser.getConfiguration() ;
    std::map<std::string, std::string> direct = parser.getDirectConfiguration() ;
    std::vector<std::string> flags = parser.getActiveFlags() ;

    ConfigTreeItem * problem = ConfigParser::readFile(file, define, true, true, flags ) ;
    problem->configure( direct ) ;

    if(printConfigTree)
    {
        problem->printTree() ;
        exit(0) ;
    }

    FeatureTree F(problem->getChild("sample")->getSample( )) ;
    if(problem->hasChildFromFullLabel("sample.sampling_number"))
        F.setSamplingFactor( F.getFeature(0), problem->getData("sample.sampling_number", 1.) ) ;
    F.setDiscretizationParameters(problem->getChild("discretization")) ;
    Vector instants = F.setSteppingParameters(problem->getChild("stepping")) ;
    InclusionFamily * allFeatures = nullptr ;
    if(problem->hasChild("inclusions"))
        allFeatures = problem->getChild("inclusions")->makeInclusionFamily( &F ) ;
    else
        allFeatures = new InclusionFamily( ) ;
    F.step() ;

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

    if(problem->hasChild("boundary_conditions"))
        std::vector<BoundaryCondition * > bc = problem->getChild("boundary_conditions")->getAllBoundaryConditions(&F) ;


    MultiTriangleWriter * trg = nullptr ;
    if(problem->hasChildFromFullLabel("export.file_name"))
    {
        std::string trgFileName = problem->getStringData("export.file_name","file_not_found") ;
        std::string headerFileName = trgFileName + "_header" ;
        trg = new MultiTriangleWriter( headerFileName, trgFileName, &F, 1.) ;
    }

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
    }

    gettimeofday ( &time1, nullptr );
    double delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
    std::cout << "problem solved in " << delta/1000000 << " seconds" << std::endl ;


    return 0 ;
}
