// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/boundarycondition.h"
#include "../features/features.h"
#include "../physics/material_laws/material_laws.h"
#include "../polynomial/vm_function_base.h"
#include "../utilities/configuration.h"
#include "../utilities/parser.h"

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>

using namespace Amie ;

int main(int argc, char *argv[])
{
	std::string file = "../examples/data/composite/test_2d_composite.ini" ;
	if(argc > 1)
		file = std::string(argv[1]) ;

	ConfigTreeItem * define = nullptr ;
	for(int i = 2 ; i < argc ; i+=2)
	{
		std::string test = std::string(argv[i]) ;
		if(test[0] == '@')
		{
			if(!define)
				define = new ConfigTreeItem(nullptr, "define") ;
			std::string testval = std::string(argv[i+1]) ;
			bool isDouble = (testval.find_first_not_of("0123456789.e-") == std::string::npos ) ;
			if(isDouble)
				ConfigTreeItem * next = new ConfigTreeItem( define, test, atof(testval.c_str()) ) ;
			else
				ConfigTreeItem * next = new ConfigTreeItem( define, test, testval ) ;
		}
	}

	ConfigTreeItem * problem = ConfigParser::readFile(file, define) ;

	FeatureTree F(problem->getChild("sample")->getSample()) ;
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
			std::vector<std::vector<Feature *> > tmp = newInclusions[i]->getInclusions( &F, dummy, inclusions ) ;
			for(size_t j = 0 ; j < tmp[0].size() ; j++)
				inclusions.push_back( dynamic_cast<Geometry *>(tmp[0][j]) ) ;
			for(size_t j = 0 ; j < tmp.size() ; j++)
			{
				std::vector<Geometry *> geom ;
				for(size_t k = 0 ; k < tmp[j].size() ; k++)
					geom.push_back( dynamic_cast<Geometry *>(tmp[j][k]) ) ;
				allFeatures.push_back( geom ) ;
			}
		}
	}

	F.step() ;
	std::vector<unsigned int> cacheIndex ;
	std::cout << "generating cache for inclusion family 0/" << allFeatures.size() ;
	cacheIndex.push_back( F.get2DMesh()->generateCache( F.getFeature(0)) ) ;
	for(size_t i = 0 ; i < allFeatures.size() ; i++)
	{
		std::cout << "\rgenerating cache for inclusion family " << i+1 << "/" << allFeatures.size() ;
		cacheIndex.push_back( F.get2DMesh()->generateCache( allFeatures[i] ) ) ;
	}
	std::cout << "... done" << std::endl ;
	for(size_t i = 0 ; i < cacheIndex.size() ; i++)
		std::cout << "inclusion family " << i << " covering surface " << F.get2DMesh()->getArea( cacheIndex[i] ) << std::endl ;

	std::vector<ConfigTreeItem *> bcItem = problem->getAllChildren("boundary_condition") ;
	std::vector<std::pair<BoundaryCondition *, LinearInterpolatedExternalMaterialLaw *> > interpolatedBC ;
	std::vector<std::pair<BoundaryCondition *, Function> > functionBC ;
	for(size_t i = 0 ; i < bcItem.size() ; i++)
	{
		BoundaryCondition * bc = bcItem[i]->getBoundaryCondition() ;
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

	for(size_t i = 1 ; i < instants.size() ; i++)
	{
		F.setDeltaTime( instants[i]-instants[i-1] ) ;

		for(size_t j = 0 ; j < interpolatedBC .size() ; j++)
			interpolatedBC[j].first->setData( interpolatedBC[j].second->get( instants[i] ) ) ;
		for(size_t j = 0 ; j < functionBC .size() ; j++)
			functionBC[j].first->setData( VirtualMachine().eval(functionBC[j].second, 0., 0., 0., instants[i] ) ) ;

		F.step() ;
		if(problem->hasChild("output"))
			problem->getChild("output")->writeOutput(&F, i, instants.size(), cacheIndex) ;
		if(problem->hasChild("export"))
			problem->getChild("export")->exportSvgTriangles(trg, &F, i, instants.size()) ;
	}
		
	return 0 ;
}
