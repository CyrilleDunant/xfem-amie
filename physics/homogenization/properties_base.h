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
	TAG_MAX_STRAIN,
	TAG_MAX_TENSILE_STRAIN,
	TAG_MAX_COMPRESSIVE_STRAIN,
	TAG_IMPOSED_STRAIN,
	TAG_STRESS,
	TAG_MAX_STRESS,
	TAG_MAX_TENSILE_STRESS,
	TAG_MAX_COMPRESSIVE_STRESS,
	TAG_RUPTURE_ENERGY,
	TAG_EXPANSION_COEFFICIENT, 
	TAG_CRACK_DENSITY, 
	TAG_CRACK_LENGTH,
	TAG_DIFFUSION_COEFFICIENT,
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
	const Tag & tag() const {return ptag ; } ;
	/* \brief checks if the property is from a specific tag */
	bool is(Tag t) const {return (t == ptag) ; } ;
	/* \brief changes the tag to NULL */
	void kill() {ptag = TAG_NULL ; } ;
	/* \brief checks if the property is NULL */
	bool isNull() const {return ptag == TAG_NULL ; } ;

	/* \brief changes the tag to a specific tag. 
	You should use the GeneralConverter instead of this method */
	void set(const Tag & t) {ptag = t ; } ;
	/* \brief changes the value of the property */
	void set(double v) {p = v ; } ;
	void set(const std::string & s) ;

	/* \brief print the property */
	void print() const ;
	void print(const std::string & indent) const ;
} ;


class Scheme ;

/**
* A material is simply a vector of properties. It is possible to fix some properties to a specific
* value with the tagset attribute. Each tag defined in the tagset will not be changed unless forced to.
*/
class Material : public std::vector<Properties>
{
protected:
	std::string name ;
	std::vector<Tag> tagset ;
	std::vector<Material> phases ;

public:
	/* \brief simple constructor */
	Material() ;
	Material(const std::string & n) ;
	/* \brief simple constructor from a pre-defined material */
	Material(PredefinedMaterial mat) ;
	/* \brief simple constructor from a core property*/
	Material(const Properties & p) ;
	/* \brief simple constructor from a vector of properties*/
	Material(const std::vector<Properties> & p) ;
	/* \brief builds a material with the mechanical properties corresponding to a CauchyGreen tensor*/
	Material(const Matrix & cauchy) ;

	/* \brief returns the list of all index for the properties of a specific tag*/
	std::vector<int> getIndex(const Tag & t) const ;
	/**
	* @param t a specific tag
	* @param i the desired position (-1 for last, 0 for first)
	* @return the index of the ith property in the vector with tag t (-1 if the tag does not exist in the material)
	*/
	int getIndex(const Tag & t, int i) const ;
	/* \brief returns the value for a specific tag (i is the index, asd defined for getIndex(Tag t, int i) */
	double val(const Tag & t, int i) const ;

	/**
	* Adds p to the material, and kills any other property with the same tag as p. If p is from a tag which is set
	* in tagset, the replacement fails.
	*/ 
	bool replace(const Properties & p) ;
	void add(const Tag & t, double val) {this->push_back(Properties(t,val)) ; } ;
	/* \brief as replace, but force the replacement even if the tag for p is set in the material */
	bool replaceForce(const Properties & p) ;
	/* \brief checks if a specific tag is set*/
	bool isSet(const Tag & t) const ;
	/* \brief set a tag*/
	bool set(const Tag & t) ;
	/* \brief set a tag, conserve only the ith tagged property (i is the index, as defined for getIndex(Tag t, int i) */
	bool set(const Tag & t, int i) ;

	bool kill(const Tag & t) ;

	/**
	* Fuse two materials together if and only if they have the same properties on the compare tags. The values on the combine tag
	* are added together
	*/
	bool combine(const Material & m, const std::vector<Tag> & compare, const Tag & combine) ;
	/**
	* adds all the properties of m to the current material
	*/
	bool merge(const Material & m) ;

	void rename(const std::string & n) {name = n ; } ;

	/* \brief tries to find missing values for all tags t if possible */
	bool findMissing(const std::vector<Tag> & t) ;
	/* \brief tries to find a value for the missing tag t, using a GeneralConverter */
	bool findMissing(const Tag & t) ;

	Material operator*(const std::string & s) ;
	Material operator+(const Material & m) ;
	Material operator-(int j) ;
	double operator() (const Tag & t) const{return val(t,-1) ; };
	double operator() (const Tag &t, double d) {replace(Properties(t,d)) ; return val(t,-1) ;} ;

	/* \brief prints all the material properties*/
	void print() const ;
	void print(const std::string & indent) const ;

	bool build(Scheme * s, bool self) ;

	void makeFraction(bool volume) ;

	int nPhases() const {return phases.size() ; } ;
	Material & child(int i) {return phases[i] ; } ;
	const Material & child(int i) const {return phases[i] ; } ;
	void add(Tag t, double val, int i) {phases[i].add(t,val) ; } ;
	void cleanComposite() {phases.clear() ; } ;

	void divide(int i, const std::vector<double> & f, bool v) ;

} ;


} ;


#endif

