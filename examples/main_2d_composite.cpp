// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
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
	F.setDiscretizationParameters(problem->getChild("discretization")) ;
	int nsteps = F.setSteppingParameters(problem->getChild("stepping")) ;
	if(problem->hasChild("inclusions"))
		problem->getChild("inclusions")->getInclusions(&F) ;

	F.step() ;

	std::vector<ConfigTreeItem *> bc = problem->getAllChildren("boundary_condition") ;
	for(size_t i = 0 ; i < bc.size() ; i++)
		F.addBoundaryCondition(bc[i]->getBoundaryCondition()) ;

	for(int i = 0 ; i < nsteps ; i++)
	{
		F.step() ;
		problem->getChild("output")->writeOutput(&F, i, nsteps) ;
		problem->getChild("export")->exportTriangles(&F, i, nsteps) ;
	}
		
	return 0 ;
}
