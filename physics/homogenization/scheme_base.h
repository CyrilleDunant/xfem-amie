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


/**
* A Scheme is an object which basically transforms materials into new entities. It is mainly used
* to derived homogenized properties of a composite material from the properties of each constituent
* of the composite, but it may be used to make other manipulations.<br>
* To implement a new scheme, you need to create a new class which inherits from Scheme. In the
* constructors of this new class, you need to make the list of all input and output properties that
* the scheme needs. You also need to define how many materials you need for the scheme (-1 for infinite).<br>
* Then you only need to overwrite the process method, writing all the formula and the algorithm of the scheme
* in this method. Don't hesitate to add new Tag if you need (in the "properties_base.h" file) or new statuses.
*/
class Scheme
{
protected:
	/* \brief describes if everything is OK in the scheme, or if an error has occured*/
	Status s ;
	/* \brief the number of phases in the composite material, -1 for infinite*/
	int p ;
	/* \brief list of input properties needed*/
	std::vector<Tag> input ;
	/* \brief list of output properties*/
	std::vector<Tag> output ;

public:
	/* \brief simple constructor*/
	Scheme(int n) ;
	/* \brief constructor, using a single input, output is the same*/
	Scheme(int n, Tag in) ;
	/* \brief constructor, using a single input and a different output*/
	Scheme(int n, Tag in, Tag out) ;
	/* \brief constructor, using a set of different input (output are the same)*/
	Scheme(int n, std::vector<Tag> & in) ;
	/* \brief constructor, using a set of different input, and another set of different output*/
	Scheme(int n, std::vector<Tag> & in, std::vector<Tag> & out) ;

	/* \brief checks if an error has not occured */
	virtual bool isOK() const {return s == STATUS_OK || s == STATUS_RESET; } ;
	/* \brief get the current status of the scheme */
	virtual Status const status() {return s ; } ;
	/* \brief reset the status */
	virtual void reset() {s = STATUS_RESET ; } ;
	/* \brief checks if an error has occured, and then reset the status*/
	virtual bool check(bool r) ;
	virtual int phases() {return p ; } ;

	/* \brief return the input list. This method is useful in cunjonction with Material::findMissing(std::vector<Tag> t) */
	virtual std::vector<Tag> inputList() {return input ; } ;
	/* \brief return the output list*/
	virtual std::vector<Tag> outputList() {return input ; } ;

	/**
	* Applies the scheme over a list of materials. There will be an error if there is not enough materials,
	* if the there is an input material which has not the needed properties, or if an error occurs during the
	* homogenization (division by 0, for example). <br>
	* This method first makes all the needed verification, then builds a Matrix data in which all the
	* input properties are sorted according to the input list (columns) and by material (rows). The matrix is then
	* processed, and a vector of properties is then build from the output of the process method.
	*/
	virtual std::vector<Properties> homogenize(const std::vector<Material> & mat) ;
	/* \brief makes the homogenization of a single material, if possible)*/
	virtual std::vector<Properties> homogenize(const Material & mat) ;
	/* \brief makes the homogenization of two materials, if possible)*/
	virtual std::vector<Properties> homogenize(const Material & m1, const Material & m2) ;

	/**
	* This method is the most important of a Scheme, since it is the application of the scheme itself.
	*/
	virtual Vector process(const Matrix & data) ;

	/* \brief helper method to check if a double is equal to zero*/
	virtual bool equalsZero(double x) ;
	/* \brief helper method to check if a double is greater than zero*/
	virtual bool lessThanZero(double x) ;
	/* \brief helper method to make a square root, verifies that it is possible, and gives an error if not*/
	virtual double simpleSquareRoot(double square) ;
	/* \brief helper method to make a division, verifies that it is possible, and gives an error if not*/
	virtual double simpleDivision(double num, double denom) ;

	/* \brief prints the current status*/
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

