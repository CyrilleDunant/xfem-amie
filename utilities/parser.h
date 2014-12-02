
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __PARSER_H_
#define __PARSER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>
#include <map>

#include "../geometry/geometry_2D.h"
#include "../geometry/geometry_3D.h"
#include "matrixops.h"
#include "configuration.h"

using namespace Amie ;

//#define DEBUG 

class Parser
{
protected:
	std::ifstream file ;
public:
	Parser(const char* f) : file(f) { }	;
	virtual ~Parser() { };
	virtual void readData() = 0;
} ;

class InclusionParser : public Parser
{
protected:
	std::vector< Circle > * circs ;
public:
InclusionParser(const char* f): Parser(f) { circs = new std::vector< Circle >() ; } 
	InclusionParser(): Parser("inc.dt") { circs = new std::vector< Circle >() ; } 
	virtual ~InclusionParser() { delete circs ;}
	
	virtual void readData() ;
	std::vector< Amie::Circle > * getData() {return circs;}
} ;

class SegmentsParser : public Parser
{
protected:
	std::vector < SegmentedLine > * segments ;

public:
SegmentsParser(const char* f): Parser(f) { segments = new std::vector< SegmentedLine >() ; }
	SegmentsParser(): Parser("ck.dt") { segments = new std::vector< SegmentedLine >() ; }
	virtual ~SegmentsParser() { delete segments ;}
	
	virtual void readData() ;
	std::vector< Amie::SegmentedLine > * getData() {return segments;}
} ;

class PointParser : public Parser
{
protected:
	std::vector< PointSet * > * pointSets ;
	
public:
	PointParser(const char* f): Parser(f) { pointSets = new std::vector<PointSet * >() ; }
	PointParser(): Parser("point.dt") { pointSets = new std::vector<PointSet * >() ; } 
	virtual ~PointParser() {  delete pointSets ; }
	
	virtual void readData() ;
	std::vector< PointSet * > * getData() { return pointSets ;}
} ;

class PointParser3D : public Parser
{
protected:
	std::vector< Point * > * pointSets ;
	std::vector< size_t > PeriodicIds;
	
public:
	PointParser3D(const char* f): Parser(f) { pointSets = new std::vector<Point * >() ; }
	PointParser3D(): Parser("cord.getT()xt") { pointSets = new std::vector<Point * >() ; } 
	virtual ~PointParser3D() { delete pointSets ; }
	
	virtual void readData() ;
	std::vector< Point * > * getData() { return pointSets ;}
	const std::vector< size_t > & getPeriodicIds() { return PeriodicIds ;}
//	std::vector< size_t > & getPeriodicIds() { return PeriodicIds ;}
} ;

class PeriodParser : public Parser
{
protected:
	std::vector<std::pair<size_t,std::valarray<int> > > * periodSets ;
	
public:
	PeriodParser(const char* f): Parser(f) { periodSets = new std::vector<std::pair<size_t,std::valarray<int> > >() ; }
	PeriodParser(): Parser("period.getT()xt") {periodSets = new std::vector<std::pair<size_t,std::valarray<int> > > ; }
	virtual ~PeriodParser() { delete periodSets ; }
	
	virtual void readData() ;
	std::vector< std::pair<size_t,std::valarray<int> > > * getData() { return periodSets ;}
} ;

class TetrahedronParser : public Parser
{
protected:
	std::vector< std::pair<std::valarray<int>, size_t> > * tets ;
	
public:
	TetrahedronParser(const char* f): Parser(f) { tets = new std::vector<std::pair<std::valarray<int>, size_t> >() ; }
	TetrahedronParser(): Parser("conec.getT()xt") { tets = new std::vector<std::pair<std::valarray<int>, size_t> >  ; } 
	virtual ~TetrahedronParser() { delete tets ; } 
	
	virtual void readData() ;
	std::vector<std::pair<std::valarray<int>, size_t> > * getData() { return tets ;}
} ;

class HexahedronParser : public Parser
{
protected:
	std::vector< std::pair<std::valarray<int>, size_t> > * hex ;
	
public:
	HexahedronParser(const char* f): Parser(f) { hex = new std::vector<std::pair<std::valarray<int>, size_t> >() ; }
	HexahedronParser(): Parser("conec.getT()xt") { hex = new std::vector<std::pair<std::valarray<int>, size_t> >  ; } 
	virtual ~HexahedronParser() { delete hex ; } 
	
	virtual void readData() ;
	std::vector<std::pair<std::valarray<int>, size_t> > * getData() { return hex ;}
} ;


class BoundaryParser : public Parser
{
protected:
	std::vector<std::valarray<int> > * tri0 ;
	std::vector<std::valarray<int> > * tri1 ;
	std::vector<std::valarray<int> > * tri2 ;
	std::vector<std::valarray<int> > * tri3 ;
	std::vector<std::valarray<int> > * tri4 ;
	std::vector<std::valarray<int> > * tri5 ;
	
public:
	BoundaryParser(const char* f): Parser(f) 
	{ 
		tri0 = new std::vector<std::valarray<int> >() ;
		tri1 = new std::vector<std::valarray<int> >() ;
		tri2 = new std::vector<std::valarray<int> >() ;
		tri3 = new std::vector<std::valarray<int> >() ;
		tri4 = new std::vector<std::valarray<int> >() ;
		tri5 = new std::vector<std::valarray<int> >() ;
	}
	BoundaryParser(): Parser("front.getT()xt") 
	{
		tri0 = new std::vector<std::valarray<int> >() ;
		tri1 = new std::vector<std::valarray<int> >() ;
		tri2 = new std::vector<std::valarray<int> >() ;
		tri3 = new std::vector<std::valarray<int> >() ;
		tri4 = new std::vector<std::valarray<int> >() ;
		tri5 = new std::vector<std::valarray<int> >() ;
	}
	virtual ~BoundaryParser()
	{
		delete tri0 ;
		delete tri1 ;
		delete tri2 ;
		delete tri3 ;
		delete tri4 ;
		delete tri5 ;
	}
	
	virtual void readData() ;
	std::vector<std::valarray<int> > * getData(size_t max_index) ;
} ;

class HexahedronBoundaryParser : public Parser
{
protected:
	std::vector<std::valarray<int> > * carre0 ;
	std::vector<std::valarray<int> > * carre1 ;
	std::vector<std::valarray<int> > * carre2 ;
	std::vector<std::valarray<int> > * carre3 ;
	std::vector<std::valarray<int> > * carre4 ;
	std::vector<std::valarray<int> > * carre5 ;
	
public:
	HexahedronBoundaryParser(const char* f): Parser(f) 
	{ 
		carre0 = new std::vector<std::valarray<int> >() ;
		carre1 = new std::vector<std::valarray<int> >() ;
		carre2 = new std::vector<std::valarray<int> >() ;
		carre3 = new std::vector<std::valarray<int> >() ;
		carre4 = new std::vector<std::valarray<int> >() ;
		carre5 = new std::vector<std::valarray<int> >() ;
	}
	HexahedronBoundaryParser(): Parser("front.getT()xt") 
	{
		carre0 = new std::vector<std::valarray<int> >() ;
		carre1 = new std::vector<std::valarray<int> >() ;
		carre2 = new std::vector<std::valarray<int> >() ;
		carre3 = new std::vector<std::valarray<int> >() ;
		carre4 = new std::vector<std::valarray<int> >() ;
		carre5 = new std::vector<std::valarray<int> >() ;
	}
	virtual ~HexahedronBoundaryParser()
	{
		delete carre0 ;
		delete carre1 ;
		delete carre2 ;
		delete carre3 ;
		delete carre4 ;
		delete carre5 ;
	}
	
	virtual void readData() ;
	std::vector<std::valarray<int> > * getData(size_t max_index) ;
} ;

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

	static ConfigTreeItem * readFile(std::string f, ConfigTreeItem * def, bool define = true) ;

} ;

#endif // __PARSER_H_
