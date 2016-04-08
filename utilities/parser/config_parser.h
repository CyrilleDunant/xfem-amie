
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

typedef enum
{
    XML_OPEN_COMMENT,
    XML_CLOSE_COMMENT,
    XML_INLINE_COMMENT,
    XML_OPEN,
    XML_CLOSE,
    XML_INLINE,
    XML_VALUE,
    XML_INVALID,
    XML_HEADER,
} XMLTokenType ;



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
	static ConfigTreeItem * readXMLFile(std::string f, ConfigTreeItem * def, bool define = true, bool bind = false, std::vector<std::string> flags = std::vector<std::string>(), std::string path = std::string()) ;
	static std::vector<BoundaryCondition *> getBoundaryConditions( std::string filename, FeatureTree * F) ;

	static Form * getBehaviour( std::string filename, Form * def = new VoidForm(), SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL ) ;

} ;


class ConfigXMLParser : public ConfigParser
{
public:
	ConfigXMLParser(std::string f, bool a = true) : ConfigParser(f,a) { }
	ConfigXMLParser(const char * f, bool a = true) : ConfigParser(f,a) { }
	ConfigXMLParser() : ConfigParser("input.xml") { }

	virtual ~ConfigXMLParser() { } ;
	virtual void readData() ;

	static std::vector<std::string> breakLine( std::string line, std::string ignore = "\"\'" ) ;
	static XMLTokenType getTokenType( std::string test ) ;
	static std::string cropToken(std::string token, XMLTokenType type ) ;
	
	ConfigTreeItem * parseXMLToken( std::string xml ) ;
	ConfigTreeItem * parseXMLAttribute( std::string attr ) ;


} ;




#endif // __PARSER_H_
