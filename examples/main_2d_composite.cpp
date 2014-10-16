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

	ConfigTreeItem * problem = ConfigParser::readFile(file) ;

	FeatureTree F(problem->getChild("sample")->getSample()) ;
	if(problem->hasChildFromFullLabel("sample.sampling_number"))
		F.setSamplingFactor( F.getFeature(0), problem->getData("sample.sampling_number", 1.) ) ;
	F.setDiscretizationParameters(problem->getChild("discretization")) ;
	Vector instants = F.setSteppingParameters(problem->getChild("stepping")) ;
	if(problem->hasChild("inclusions"))
	{
		std::vector<Geometry *> inclusions ;
		std::vector<Feature *> dummy ;
		std::vector<ConfigTreeItem *> newInclusions = problem->getAllChildren("inclusions") ;
		for(size_t i = 0 ; i < newInclusions.size() ; i++)
		{
			std::vector<Feature *> tmp = newInclusions[i]->getInclusions( &F, dummy, inclusions ) ;
			for(size_t j = 0 ; j < tmp.size() ; j++)
				inclusions.push_back( dynamic_cast<Geometry *>(tmp[j]) ) ;
		}
	}

	F.step() ;

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
		problem->getChild("output")->writeOutput(&F, i, instants.size()) ;
		problem->getChild("export")->exportSvgTriangles(trg, &F, i, instants.size()) ;
	}
		
	return 0 ;
}
