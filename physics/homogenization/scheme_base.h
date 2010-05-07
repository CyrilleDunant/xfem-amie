//
// C++ Interface: mechanical homogenization
//
// Description: 
//
//
// Author:  Alain Giorla, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SCHEME_BASE_H
#define SCHEME_BASE_H


#include "properties_base.h"

namespace Mu
{

typedef enum
{
	STATUS_RESET,
	STATUS_OK,
	STATUS_MATERIAL_NOT_FOUND,
	STATUS_PROPERTIES_NOT_FOUND,
	STATUS_BAD_HOMOGENIZATION,
} Status ;


/* Default homogenization scheme. The actual schemes are inherited from this class*/
class Scheme
{
protected:
	Status s ;
	size_t p ;
	std::vector<Tag> input ;
	std::vector<Tag> output ;

public:
	Scheme(size_t n) ;
	/* \brief constructor, using a single input, output is the same*/
	Scheme(size_t n, Tag in) ;
	/* \brief constructor, using a single input and a different output*/
	Scheme(size_t n, Tag in, Tag out) ;
	/* \brief constructor, using a set of different input (output are the same)*/
	Scheme(size_t n, std::vector<Tag> & in) ;
	/* \brief constructor, using a set of different input, and another set of different output*/
	Scheme(size_t n, std::vector<Tag> & in, std::vector<Tag> & out) ;

	virtual bool isOK() const {return s == STATUS_OK ; } ;
	virtual Status const status() {return s ; } ;
	virtual void reset() {s = STATUS_RESET ; } ;
	virtual bool check(bool r) ;

	virtual std::vector<Tag> inputList() {return input ; } ;
	virtual std::vector<Tag> outputList() {return input ; } ;

	virtual std::vector<Properties> homogenize(const std::vector<Material> & mat) ;
	virtual std::vector<Properties> homogenize(const Material & mat) ;
	virtual Vector process(const Matrix & data) ;

	virtual bool equalsZero(double x) ;
	virtual bool lessThanZero(double x) ;
	virtual double simpleSquareRoot(double square) ;
	virtual double simpleDivision(double num, double denom) ;

	virtual void print() ;

} ;




class MeanScheme : public Scheme
{
protected:
	bool parallel ;

public:
	MeanScheme(bool volume, bool parallel, Tag t) ;
	MeanScheme(bool volume, bool parallel, std::vector<Tag> t) ;

	virtual Vector process(const Matrix & data) ;
	virtual Vector processParallel(const Matrix & data) ;

	
} ;




} ;


#endif // SCHEME

