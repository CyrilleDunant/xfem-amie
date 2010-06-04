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

typedef enum
{
	CIRCLE_INCLUSION,
	SPHERE_INCLUSION,
	ELLIPSE_INCLUSION
} TypeInclusion;

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

	virtual std::vector <Inclusion *> operator()(double , double, int inclusionNumber = 8000, double itzSize = 15e-6, bool verbose = true);

	virtual std::vector <EllipsoidalInclusion *> operator()(bool, Point *, double , double, double rfactor = 0.8, int inclusionNumber = 8000, double itzSize = 15e-6, bool verbose = true);
} ;

class GranuloBolome : public Granulo
{
private:
	TypeGranulo type ;
	double masseInitiale;
	double densite;
public:
	GranuloBolome(double , double , TypeGranulo t) ;
	virtual std::vector <Inclusion *> operator()(double , double, double inclusionNumber = 8000, double itzSize = 15e-6, bool verbose = true);
	virtual std::vector <Inclusion3D *> operator()(bool, double rayonGranulatMax, double pourcentMasseMin, double inclusionNumber = 8000, double itzSize = 15e-6, bool verbose = true);

	virtual ~GranuloBolome() { } ;
} ;

class GranuloFromFile
{
private:
	std::string filename;
	std::vector<std::string> fields;
	std::vector<double> values ;

public:
	GranuloFromFile(std::string fname, std::vector<std::string> columns) ;
	int numberOfField() {return this->fields.size() ; } ;
	bool verifyField(std::vector<std::string> columns) ;
	int getFieldNumber(std::string column) ;
	std::vector<double> getFieldValues(std::string column) ;
	std::vector<double> getFieldValues(int) ;
	std::vector<Feature *> getFeatures(TypeInclusion type, int ninc) ;
	std::vector<Inclusion3D *> getInclusion3D(int ninc) ;
} ;



}

