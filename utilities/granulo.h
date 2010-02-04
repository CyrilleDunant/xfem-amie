//
// C++ Interface: granulo
//
// Description: 
//
//
// Author:  <>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <iostream>
#include <vector>
#include "../features/inclusion.h"
#include "../features/inclusion3d.h"
#include "../geometry/geometry_base.h"


namespace Mu
{

typedef enum
{
	BOLOME_A,
	BOLOME_B,
	BOLOME_C,
	BOLOME_D
} TypeGranulo;

class Granulo
{	
private:
	double c;
	double n;
	double masseInitiale;
	double densite;
public:

	Granulo(double, double, double, double);
	Granulo() ;
	virtual ~Granulo() { } ;

	virtual std::vector <Inclusion *> operator()(double , double, int inclusionNumber = 8000, double itzSize = 15e-6);

	virtual std::vector <EllipsoidalInclusion *> operator()(bool, Point *, double , double, double rfactor = 0.8, int inclusionNumber = 8000, double itzSize = 15e-6);
} ;

class GranuloBolome : public Granulo
{
private:
	TypeGranulo type ;
	double masseInitiale;
	double densite;
public:
	GranuloBolome(double , double , TypeGranulo t) ;
	virtual std::vector <Inclusion *> operator()(double , double, int inclusionNumber = 8000, double itzSize = 15e-6);
	virtual std::vector <Inclusion3D *> operator()(bool, double , double, int inclusionNumber = 8000, double itzSize = 15e-6);

	virtual ~GranuloBolome() { } ;
} ;

class GranuloFromFile
{
private:
	std::string filename;
	std::vector<double> mass ;
	std::vector<double> size ;

public:
	GranuloFromFile(std::string fname, std::vector<double>, double, double) ;
	void resize(double) ;
	std::vector<Inclusion *> getCircleInclusion(double,int,double) ;
	std::vector<EllipsoidalInclusion *> getEllipsoidalInclusion(double,int,double) ;
} ;



}

