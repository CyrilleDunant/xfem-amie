
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __CONFIG_PARSER_H_
#define __CONFIG_PARSER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>
#include <map>

#include "../../geometry/geometry_2D.h"
#include "../../geometry/geometry_3D.h"
#include "../matrixops.h"
#include "../configuration.h"
#include "parser.h"

using namespace Amie ;

//#define DEBUG 



class ConfigParser : public Parser
{
protected:
	ConfigTreeItem * trunk ;
	std::string filename ;
	bool authorizeIncludes ;

	static int getIndentLevel( std::string test ) ;

public:
	ConfigParser(std::string f, bool a = true): Parser(f.c_str()), filename(f), authorizeIncludes(a) { trunk = new ConfigTreeItem() ; } 
	ConfigParser(const char* f, bool a = true): Parser(f), filename(f), authorizeIncludes(a) { trunk = new ConfigTreeItem() ; } 
	ConfigParser(): Parser("input.ini"), filename("input.ini") { trunk = new ConfigTreeItem() ; } 
	virtual ~ConfigParser() { } ;
	virtual void readData() ;
	ConfigTreeItem * getData() { return trunk; }

	static ConfigTreeItem * readFile(std::string f, ConfigTreeItem * def, bool define = true, bool bind = false, std::vector<std::string> flags = std::vector<std::string>(), std::string path = std::string()) ;
	static std::vector<BoundaryCondition *> getBoundaryConditions( std::string filename, FeatureTree * F) ;

	static Form * getBehaviour( std::string filename, Form * def = new VoidForm(), SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL ) ;

} ;



#endif // __PARSER_H_
