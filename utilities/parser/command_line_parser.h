
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __COMMAND_LINE_PARSER_H_
#define __COMMAND_LINE_PARSER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>
#include <map>

#include "../configuration.h"

using namespace Amie ;

//#define DEBUG 

struct CommandLineArgument
{
	std::string name ;
	std::string help ;
	std::string str ;
	double val ;

	CommandLineArgument(std::string n, std::string h = std::string(), std::string defstr = std::string(), double v = 0.) : name(n), help(h), str(defstr), val(v) { } ;
} ;

class CommandLineParser
{
protected:
	std::map< std::string, bool > flags ;
	std::map< std::string, double > values ;
	std::map< std::string, std::string > strings ;
	std::map< std::string, std::string> help ;
	std::vector< CommandLineArgument > arguments ;
	std::map<std::string, std::string> aliases ;
	std::vector<std::string> parsed ;
	std::string command ;
	std::string commandLine ;
	std::string description ;

	bool commandLineConfiguration ;
	bool forceUnrecognizedFlags ;
	ConfigTreeItem * config ;
	ConfigTreeItem * input ;
	std::map<std::string, std::string> directConfig ;

public:
	CommandLineParser(std::string d = std::string(), bool c = false, bool f = false) ;

	void addAlias(std::string alias, std::string complete) { aliases[alias] = complete ; }
	void addFlag( std::string f, std::string h = std::string(), std::string a = std::string() ) { flags[f] = false ; help[f] = h ; if(a.length() > 0) { addAlias(a,f) ; } }
	void addValue( std::string f, double val, std::string h = std::string(), std::string a = std::string() ) { values[f] = val ; help[f] = h ; if(a.length() > 0) { addAlias(a,f) ; } }
	void addString( std::string f, std::string str, std::string h = std::string(), std::string a = std::string() ) { strings[f] = str ; help[f] = h ; if(a.length() > 0) { addAlias(a,f) ; } }
	void addArgument( std::string f, std::string defstr, std::string h = std::string()) { arguments.push_back( CommandLineArgument(f, h, defstr, 0. ) ) ; }
	void addArgument( std::string f, double v, std::string h = std::string()) { arguments.push_back( CommandLineArgument(f, h, std::string(), v ) ) ; }
	std::string getCompleteString(std::string alias) ;
	std::string getAlias(std::string complete) ;

	bool getFlag( std::string f ) ;
	double getValue( std::string f ) ;
	std::string getString( std::string f ) ;
	std::string getStringArgument( std::string arg ) ;
	std::string getStringArgument( size_t i ) ;
	double getNumeralArgument( std::string arg ) ;
	double getNumeralArgument( size_t i ) ;
	Form * getBehaviour( std::string arg, Form * b, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;

	ConfigTreeItem * parseCommandLine( int argc, char *argv[], std::vector<std::string> sargs = std::vector<std::string>() ) ;
	ConfigTreeItem * parseConfigFile( std::string file, bool priority = false ) ;

	std::vector<std::string> getActiveFlags() ;
	std::string getCommandLine() const { return commandLine ; }

	ConfigTreeItem * getLocalConfiguration() { return config ; }
	std::map<std::string, std::string> getDirectConfiguration() { return directConfig ; }

	void setNumThreads(int n) ;
	void sendEmail( std::string subject, std::string body, std::string attachment = std::string() ) ;

	void printStatus() ;
	void printHelp() ;
	void printVersion() ;
	void printFormatedHelp( std::string arg, int max, int maxalias, std::string help, std::string lead, bool printAlias) ;

	void setFeatureTree( FeatureTree * f) ;
	static void setFeatureTree( FeatureTree * f, int argc, char *argv[], std::string str = std::string("AMIE") ) ;
	void disableFeatureTreeArguments() ;

} ;




#endif // __COMMAND_LINE_PARSER_H_
