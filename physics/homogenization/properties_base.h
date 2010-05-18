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

namespace Mu
{


typedef enum
{
	TAG_NULL, 
	TAG_UNIVERSAL, 
	TAG_VOLUME,	
	TAG_VOLUME_FRACTION, 
	TAG_VOLUME_TOTAL,
	TAG_MASS,
	TAG_MASS_FRACTION,
	TAG_MASS_TOTAL,
	TAG_DENSITY,
	TAG_YOUNG_MODULUS,
	TAG_POISSON_RATIO,
	TAG_BULK_MODULUS,	
	TAG_SHEAR_MODULUS,
	TAG_LAME_COEFFICIENT,
	TAG_STRAIN,
	TAG_STRESS,
	TAG_EXPANSION_COEFFICIENT, 
	TAG_CRACK_DENSITY, 
	TAG_CRACK_LENGTH,
	TAG_ELLIPSE_A,
	TAG_ELLIPSE_B,
	TAG_AREA,
	TAG_PERIMETER,
	TAG_ELLIPSE_FIRST_COMPLETE_INTEGRAL,
	TAG_ELLIPSE_SECOND_COMPLETE_INTEGRAL,
	TAG_CIRCLE_RADIUS,
} Tag ;

typedef enum
{
	MAT_DUMMY,
	MAT_AGGREGATE,
	MAT_CEMENT,
} PredefinedMaterial ;

/**
* A Properties item is simply a double tagged with a property type.
*/
class Properties
{
protected:
	Tag ptag ;
	double p ;

public:
	/* \brief default constructor for a NULL properties item*/
	Properties() ;
	/* \brief constructor for an UNIVERSAL properties item */
	Properties(double v) ;
	/* \brief constructor for an specific properties item */
	Properties(Tag p, double v) ;
	/* \brief copy constructor */
	Properties(const Properties & p) ;

	/* \brief returns the value of the property */
	double val() const {return p ; } ;
	/* \brief returns the tag of the property */
	Tag tag() const {return ptag ; } ;
	/* \brief checks if the property is from a specific tag */
	bool is(Tag t) const {return ptag == t ; } ;
	/* \brief changes the tag to NULL */
	void kill() {ptag = TAG_NULL ; } ;
	/* \brief checks if the property is NULL */
	bool isNull() const {return ptag == TAG_NULL ; } ;

	/* \brief changes the tag to a specific tag. 
	You should use the GeneralConverter instead of this method */
	void set(Tag t) {ptag = t ; } ;
	/* \brief changes the value of the property */
	void set(double v) {p = v ; } ;

	/* \brief print the property */
	void print() ;
} ;


/**
* A material is simply a vector of properties. It is possible to fix some properties to a specific
* value with the tagset attribute. Each tag defined in the tagset will not be changed unless forced to.
*/
class Material : public std::vector<Properties>
{
protected:
	std::vector<Tag> tagset ;

public:
	/* \brief simple constructor */
	Material() ;
	/* \brief simple constructor from a pre-defined material */
	Material(PredefinedMaterial mat) ;
	/* \brief simple constructor from a core property*/
	Material(const Properties & p) ;
	/* \brief simple constructor from a vector of properties*/
	Material(const std::vector<Properties> & p) ;
	/* \brief builds a material with the mechanical properties corresponding to a CauchyGreen tensor*/
	Material(const Matrix & cauchy) ;

	/* \brief returns the list of all index for the properties of a specific tag*/
	std::vector<int> getIndex(Tag t) const ;
	/**
	* @param t a specific tag
	* @param i the desired position (-1 for last, 0 for first)
	* @return the index of the ith property in the vector with tag t (-1 if the tag does not exist in the material)
	*/
	int getIndex(Tag t, int i) const ;
	/* \brief returns the value for a specific tag (i is the index, asd defined for getIndex(Tag t, int i) */
	double val(Tag t, int i) const ;

	/**
	* Adds p to the material, and kills any other property with the same tag as p. If p is from a tag which is set
	* in tagset, the replacement fails.
	*/ 
	bool replace(Properties p) ;
	/* \brief as replace, but force the replacement even if the tag for p is set in the material */
	bool replaceForce(Properties p) ;
	/* \brief checks if a specific tag is set*/
	bool isSet(Tag t) const ;
	/* \brief set a tag*/
	bool set(Tag t) ;
	/* \brief set a tag, conserve only the ith tagged property (i is the index, as defined for getIndex(Tag t, int i) */
	bool set(Tag t, int i) ;

	/**
	* Fuse two materials together if and only if they have the same properties on the compare tags. The values on the combine tag
	* are added together
	*/
	bool combine(Material m, std::vector<Tag> compare, Tag combine) ;
	/**
	* adds all the properties of m to the current material
	*/
	bool merge(Material m) ;

	/* \brief tries to find missing values for all tags t if possible */
	bool findMissing(std::vector<Tag> t) ;
	/* \brief tries to find a value for the missing tag t, using a GeneralConverter */
	bool findMissing(Tag t) ;

	/* \brief prints all the material properties*/
	void print() ;
} ;


} ;


#endif

