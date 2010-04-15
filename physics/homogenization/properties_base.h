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

#ifndef PROPERTIES_BASE_H
#define PROPERTIES_BASE_H

#include "../../utilities/matrixops.h"
#include "../../geometry/geometry_base.h"


namespace Mu
{


typedef enum
{
	VOID_PROP,
	ABSTRACT,
	FRACTION,
	HOOKE,
	BULK_SHEAR,
	EXPANSION,
} PropertiesType ;

/* \brief returns the standard number of value for a specific type */
size_t standardNVal(Mu::PropertiesType p) ;
/* \brief returns if the conversion is possible between p1 and p2*/
bool conversionPossible(Mu::PropertiesType p1, Mu::PropertiesType p2) ;

/* \brief Properties item are a vector (values) specified by a type
* The number of values in the vector is given by the type.
* VOID_PROP: void properties item (no values inside)
* ABSTRACT: generic properties item, any number of values
* FRACTION: used for mass or volume fraction of phases in the composite material (1 value)
* HOOKE: Young Modulus and Poisson Ratio (2 values). Can be converted to BULK_SHEAR
* BULK_SHEAR: Bulk Modulus and Shear Modulus (2 values). Can be converted to HOOKE
*/
class Properties
{
protected:
	PropertiesType pType ;
	size_t nVal ;
	Vector values ;

public:
	/* \brief default constructor for a VOID_PROP properties item*/
	Properties() ;
	/* \brief constructor for an ABSTRACT properties item */
	Properties(double v) ;
	/* \brief constructor for an ABSTRACT properties item */
	Properties(const std::pair<double, double> & v) ;
	/* \brief constructor for an ABSTRACT properties item */
	Properties(const std::vector<double> & v) ;
	/* \brief constructor for an ABSTRACT properties item */
	Properties(const Vector & v) ;
	/* \brief constructor for an specific properties item */
	Properties(PropertiesType p, double v) ;
	/* \brief constructor for an specific properties item */
	Properties(PropertiesType p, const std::pair<double, double> & v) ;
	/* \brief constructor for an specific properties item */
	Properties(PropertiesType p, const std::vector<double> & v) ;
	/* \brief constructor for an specific properties item */
	Properties(PropertiesType p, const Vector & v) ;
	/* \brief constructor for a specific properties item from a matrix. Can convert Cauchy-Green matrixes into E-nu or k-mu variables */
	Properties(PropertiesType p, const Matrix & m) ;
	/* \brief copy constructor */
	Properties(const Properties & p) ;

	/* \brief returns the actual number of values in the properties item*/
	size_t getNVal() const {return nVal ; } ;
	/* \brief returns the ith value*/
	double getValue(size_t i) const {return values[i] ; } ;
	/* \brief returns all values*/
	const Vector & getValues() const {return values ; } ;
	/* \brief returns the type of the item*/
	void setValue(size_t i, double d) {values[i] = d ; } ;
	PropertiesType getPropertiesType() const {return pType ; } ;
	std::pair<bool,Matrix> getCauchyGreen(SpaceDimensionality dim) const ;

	/* \brief converts to another type. Returns false and a VOID_PROP item if conversion is impossible*/
	std::pair<bool, Properties> convert(PropertiesType p_out) ;

	void print() ;

} ;


/* \brief vector of Properties */
class Material : public std::vector<Properties>
{
public:
	Material() ;
	Material(const Properties & p) ;
	Material(const std::vector<Properties> & p) ;

	/*  */
	size_t getFirstIndex(PropertiesType p) const ;
	/* return the position of the HOOKE properties (creates it if needed and possible) */
	size_t getHooke() ;
	/* return the position of the BULK_SHEAR properties (creates it if needed and possible) */
	size_t getBulkShear() ;

} ;


} ;


#endif

