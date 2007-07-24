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

	virtual std::vector <Inclusion *> operator()(double ,double);
} ;

class GranuloBolome : public Granulo
{
private:
	TypeGranulo type ;
	double masseInitiale;
	double densite;
public:
	GranuloBolome(double , double , TypeGranulo t) ;
	virtual std::vector <Inclusion *> operator()(double , double);
	
	virtual ~GranuloBolome() { } ;
} ;

}

