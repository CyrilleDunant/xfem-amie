
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __FUNCTION_PARSER_H_
#define __FUNCTION_PARSER_H_

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

using namespace Amie ;

//#define DEBUG 



struct FunctionParser
{
	std::vector<FunctionParser *> tokens ;
	std::vector<int> roots ;
	bool isLinking ;

	FunctionParser(std::vector<std::string> data, std::map<std::string, char> coordinates = std::map<std::string, char>()) ;
	FunctionParser(std::string f, std::map<std::string, char> coordinates = std::map<std::string, char>()) ;
	FunctionParser() : isLinking(false) { } ;

	void renewExpression( std::vector<std::string> f, std::map<std::string, char> coordinates  = std::map<std::string, char>()) ;

	std::string translateCoordinate( std::string test, std::map<std::string, char> coordinates) ;

	virtual bool isFinalToken() const { return false ; }

	virtual void print(bool end = true) const ;
	virtual void printRoots() const ;

	Function getFunction( size_t i ) const ;
	Function getFunction() const ;

	void link() ;
	void linkLeftAndRightToken( TokenOperationType op) ;
//	void linkBinaryToken( TokenOperationType op) ;

	int getLeftTokenIndex( size_t i ) const ;
	int getRightTokenIndex( size_t i ) const ;
	Function getLeftFunction( size_t i ) const ;
	Function getRightFunction( size_t i ) const ;

	static Function getFunction( std::string f,  std::map<std::string, char> coordinates = std::map<std::string, char>() ) ;
	static std::vector<std::string> breakString( std::string f ) ;
} ;

struct FunctionParserToken : public FunctionParser
{
	std::string data ;

	FunctionParserToken(std::string d) : data(d) { } ;

	virtual bool isFinalToken() const { return true ; }

	virtual void print(bool end = true ) const ;

	TokenOperationType toFunctionToken() const ;

} ;





#endif // __PARSER_H_
