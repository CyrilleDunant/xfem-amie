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

    CommandLineParser parser("Run a 2D AMIE simulation using the configuration parameters found in a *.ini file", true) ;
    parser.addFlag( "--print-microstructure", false , "print only the mesh and values of mechanical properties") ;
    parser.addArgument( "file_name", "../examples/data/composite/test_2d_composite.ini", "relative path to *.ini file to run") ;
    parser.parseCommandLine( argc, argv ) ;
    bool printMicrostructureOnly = parser.getFlag("--print-microstructure") ;
    std::string file = parser.getStringArgument(0) ;

    ConfigTreeItem * define = parser.getConfiguration() ;
    std::map<std::string, std::string> direct = parser.getDirectConfiguration() ;

    ConfigTreeItem * problem = ConfigParser::readFile(file, define) ;
    problem->configure( direct ) ;

    std::vector<ExternalMaterialLaw *> common ;
    if(problem->hasChild("common_fields"))
    {
        std::vector<ConfigTreeItem *> mat = problem->getChild("common_fields")->getAllChildren("material_law") ;
        for(size_t i = 0 ; i < mat.size() ; i++)
        {
            ExternalMaterialLaw * law = mat[i]->getExternalMaterialLaw() ;
            if(law != nullptr)
               common.push_back(law) ;
        }
    }

    FeatureTree F(problem->getChild("sample")->getSample( common )) ;
    if(problem->hasChildFromFullLabel("sample.sampling_number"))
        F.setSamplingFactor( F.getFeature(0), problem->getData("sample.sampling_number", 1.) ) ;
    F.setDiscretizationParameters(problem->getChild("discretization")) ;
    Vector instants = F.setSteppingParameters(problem->getChild("stepping")) ;
    std::vector<std::vector<Geometry *> > allFeatures ;
    if(problem->hasChild("inclusions"))
    {
        std::vector<Geometry *> inclusions ;
        std::vector<Feature *> dummy ;
        std::vector<ConfigTreeItem *> newInclusions = problem->getAllChildren("inclusions") ;
        for(size_t i = 0 ; i < newInclusions.size() ; i++)
        {
            std::vector<std::vector<Feature *> > tmp = newInclusions[i]->getInclusions( &F, dummy, inclusions, common ) ;
            for(size_t j = 0 ; j < tmp[0].size() ; j++)
                inclusions.push_back( dynamic_cast<Geometry *>(tmp[0][j]) ) ;
            for(size_t j = 0 ; j < tmp.size() ; j++)
            {
                std::vector<Geometry *> geom ;
                for(size_t k = 0 ; k < tmp[j].size() ; k++)
                {
                    if( F.getFeature(0)->in(tmp[j][k]->getCenter()) || F.getFeature(0)->intersects(tmp[j][k]))
                        geom.push_back( dynamic_cast<Geometry *>(tmp[j][k]) ) ;
                }
                allFeatures.push_back( geom ) ;
            }
        }
    }

    F.step() ;

    std::vector<unsigned int> cacheIndex ;
    std::vector<unsigned int> aggCacheIndex ;
    std::cout << "generating cache for inclusion family 0/" << allFeatures.size() ;
    cacheIndex.push_back( 0 ) ;
    for(size_t i = 0 ; i < allFeatures.size() ; i++)
    {
        std::cout << "\rgenerating cache for inclusion family " << i+1 << "/" << allFeatures.size() ;
        cacheIndex.push_back( F.get2DMesh()->generateCache( allFeatures[i] ) ) ;
        aggCacheIndex.push_back( cacheIndex[cacheIndex.size()-1] ) ;
    }
    std::cout << "generating cache for inclusion family 0/" << allFeatures.size() ;
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

    std::vector<ConfigTreeItem *> bcItem = problem->getAllChildren("boundary_condition") ;
    std::vector<std::pair<BoundaryCondition *, LinearInterpolatedExternalMaterialLaw *> > interpolatedBC ;
    std::vector<std::pair<BoundaryCondition *, Function> > functionBC ;
    for(size_t i = 0 ; i < bcItem.size() ; i++)
    {
        BoundaryCondition * bc = bcItem[i]->getBoundaryCondition(&F) ;
        if(bcItem[i]->hasChild("time_evolution"))
        {
            if(bcItem[i]->hasChildFromFullLabel("time_evolution.file_name"))
            {
                LinearInterpolatedExternalMaterialLaw * interpolation = new LinearInterpolatedExternalMaterialLaw(std::make_pair("t","value"), bcItem[i]->getStringData("time_evolution.file_name", "file_not_found")) ;
                interpolatedBC.push_back(std::make_pair(bc, interpolation)) ;
            }
            if(bcItem[i]->hasChildFromFullLabel("time_evolution.function"))
            {
                Function f = bcItem[i]->getChildFromFullLabel("time_evolution.function")->getFunction() ;
                functionBC.push_back(std::make_pair(bc, f)) ;
            }
            if(bcItem[i]->hasChildFromFullLabel("time_evolution.rate"))
            {
                Function f = "t" ;
                f *= bcItem[i]->getData("time_evolution.rate", 0.) ;
                functionBC.push_back(std::make_pair(bc, f)) ;
            }
            if(bcItem[i]->hasChildFromFullLabel("time_evolution.instants") && bcItem[i]->hasChildFromFullLabel("time_evolution.values"))
            {
                Vector t_ = ConfigTreeItem::readLineAsVector(bcItem[i]->getStringData("time_evolution.instants","0,1")) ;
                Vector v_ = ConfigTreeItem::readLineAsVector(bcItem[i]->getStringData("time_evolution.values","0,0")) ;
                LinearInterpolatedExternalMaterialLaw * interpolation = new LinearInterpolatedExternalMaterialLaw(std::make_pair("t","value"), std::make_pair(t_,v_) ) ;
                interpolatedBC.push_back(std::make_pair(bc, interpolation)) ;
            }
        }
        F.addBoundaryCondition(bc) ;
    }

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

        for(size_t j = 0 ; j < interpolatedBC .size() ; j++)
            interpolatedBC[j].first->setData( interpolatedBC[j].second->get( instants[i] ) ) ;
        for(size_t j = 0 ; j < functionBC .size() ; j++)
            functionBC[j].first->setData( VirtualMachine().eval(functionBC[j].second, 0., 0., 0., instants[i] ) ) ;

        next = F.step() ;
//		F.getAssembly()->print() ;
//		exit(0) ;
        if(problem->hasChild("output"))
            problem->getChild("output")->writeOutput(&F, i, instants.size(), cacheIndex) ;
        if(problem->hasChild("export"))
            problem->getChild("export")->exportSvgTriangles(trg, &F, i, instants.size()) ;
    }

    gettimeofday ( &time1, nullptr );
    double delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
    std::cout << "problem solved in " << delta/1000000 << " seconds" << std::endl ;


    return 0 ;
}
